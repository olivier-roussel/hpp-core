// Copyright (c) 2016, Joseph Mirabel
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
//
// This file is part of hpp-core.
// hpp-core is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-core is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-core. If not, see <http://www.gnu.org/licenses/>.

#define HPP_DEBUG 1
#define HPP_ENABLE_BENCHMARK 1

#include "hpp/core/path-optimization/recursive-hermite.hh"

#include <hpp/util/timer.hh>
#include <hpp/pinocchio/util.hh>

#include <hpp/core/path-vector.hh>
#include <hpp/core/path/hermite.hh>
#include <hpp/core/path/composite-bezier.hh>
#include <hpp/core/path-validation.hh>
#include <hpp/core/path-validation-report.hh>
#include <hpp/core/interpolated-path.hh>
#include <hpp/core/config-projector.hh>
#include <hpp/core/steering-method/hermite.hh>
#include <hpp/core/straight-path.hh>

#include <limits>

namespace hpp {
  namespace core {
    HPP_DEFINE_TIMECOUNTER(validation);

    namespace pathOptimization {
      typedef std::vector<path::HermitePtr_t> HermitePaths_t;
      typedef std::vector<HermitePaths_t> HermitePathss_t;
      const value_type infty = std::numeric_limits<value_type>::infinity();

      PathVectorPtr_t cleanInput (const PathVectorPtr_t& input)
      {
        PathVectorPtr_t flat = PathVector::create
          (input->outputSize(), input->outputDerivativeSize());
        input->flatten(flat);
        // Remove zero length path
        PathVectorPtr_t clean = PathVector::create
          (input->outputSize(), input->outputDerivativeSize());
        for (std::size_t i = 0; i < flat->numberPaths(); ++i) {
          PathPtr_t p = flat->pathAtRank (i);
          if (p->length() > 0) clean->appendPath (p);
        }
        return clean;
      }

      ConfigProjectorPtr_t getConfigProj (const PathPtr_t& p)
      {
        const ConstraintSetPtr_t& c = p->constraints();
        if (c) return c->configProjector();
        else   return ConfigProjectorPtr_t();
      }

      bool configAndVelocityAt (PathPtr_t p, const value_type& t, Configuration_t& q, vector_t& v)
      {
        if (!p) return false;
        bool success = (*p) (q, t);
        p->derivative(v, t, 1);

        ConfigProjectorPtr_t proj = getConfigProj (p);
        if (proj) proj->projectVectorOnKernel (q, v, v);
        return success;
      }

      bool velocityAt (PathPtr_t p, const value_type& t, Configuration_t& q, vector_t& v)
      {
        if (!p) return false;
        p->derivative(v, t, 1);

        ConfigProjectorPtr_t proj = getConfigProj (p);
        if (proj) {
          (*p) (q, t);
          proj->projectVectorOnKernel (q, v, v);
          return true;
        }
        return false;
      }

      /// \tparam beginning whether to use initial or final velocity.
      /// \return true if the velocity was projected.
      template <bool beginning>
      bool velocity (PathPtr_t p, Configuration_t& q, vector_t& v)
      {
        if (!p) return false;
        const interval_t& tr = p->timeRange();
        value_type t = (beginning ? tr.first : tr.second);
        return velocityAt (p, t, q, v);
      }

      void velocity (PathPtr_t before, PathPtr_t after, Configuration_t& q,
          vector_t& vb, vector_t& va, vector_t& v)
      {
        assert (before || after);
        bool beforeProjected = velocity<false> (before, q, vb),
             afterProjected  = velocity<true > (after , q, va);

        if (before && after) {
          if      ( beforeProjected &&  afterProjected) v = (vb+va)/2;
          else if (!beforeProjected &&  afterProjected) v =     va   ;
          else if ( beforeProjected && !afterProjected) v =  vb      ;
          else if (!beforeProjected && !afterProjected) v = (vb+va)/2;
        }
        else if (before)     v = vb;
        else if (after )     v = va;
      }

      InterpolatedPathPtr_t project (
          const Problem& problem,
          const PathPtr_t& path,
          const value_type& distThr)
      {
        InterpolatedPathPtr_t res = InterpolatedPath::create (
            problem.robot(),
            path->initial(), path->end(),
            path->length()
            );
           // , path->constraints());

        const Distance& distance = *problem.distance();
        const PathValidationPtr_t& validation = problem.pathValidation();

        const interval_t& tr = path->timeRange();
        Configuration_t q1 (path->outputSize()),
                        q2 (path->outputSize());
        value_type t1 = tr.first, dt = (tr.second-tr.first)/2;

        q1 = path->initial();
        while (t1 + dt < tr.second) {
          bool success = (*path) (q2, t1+dt);
          if (!success) {
            hppDout (info, "Could not project configuration at param " << t1+dt);
          }
          value_type dist = distance(q1, q2);
          if (!success || dist > distThr) {
            dt /= 2;
            continue;
          }
          // Add path validation
          StraightPathPtr_t p = StraightPath::create (problem.robot(), q1, q2, dist);
          PathPtr_t validPart;
	  PathValidationReportPtr_t report;
          if (!problem.pathValidation ()->validate (p, false, validPart, report))
          {
            dt /= 2;
            continue;
          }
          t1 += dt;
          q1 = q2;
          res->insert(t1, q1);
          dt *= 2;
        }
        return res;
      }

      PathVectorPtr_t project2 (
          const Problem& problem,
          const steeringMethod::HermitePtr_t& sm,
          const PathPtr_t& path,
          const value_type& distThr)
      {
        PathVectorPtr_t res = PathVector::create
          (path->outputSize (), path->outputDerivativeSize ());

        const PathValidationPtr_t& validation = problem.pathValidation();

        const interval_t& tr = path->timeRange();
        Configuration_t q1 (path->outputSize()),
                        q2 (path->outputSize()),
                        qt (path->outputSize());
        vector_t v1 (path->outputDerivativeSize()),
                 v2 (path->outputDerivativeSize());
        value_type t1 = tr.first, dt = (tr.second-tr.first)/2;

        q1 = path->initial();
        velocityAt (path, t1, qt, v1);
        while (true) {
          if (!(*path) (q2, t1+dt)) {
            hppDout (info, "Could not project configuration at param " << t1+dt);
            dt /= 2;
            continue;
          }
          velocityAt (path, t1+dt, qt, v2);
          path::HermitePtr_t hp = sm->steer(q1, q2, v1, v2, dt);
          value_type dist = hp->hermiteLength();
          if (dist > distThr) {
            dt /= 2;
            continue;
          }
          PathPtr_t validPart;
	  PathValidationReportPtr_t report;
          if (!problem.pathValidation ()->validate (hp, false, validPart, report))
          {
            dt /= 2;
            continue;
          }
          t1 += dt;
          q1 = q2;
          v1 = v2;
          res->appendPath(hp);
          dt *= 2;
          if (t1 >= tr.second) break;
          if (t1 + dt >= tr.second) dt = tr.second - t1;
        }
        return res;
      }

      typedef path::CompositeBezier<3> BezierCurves_t;
      typedef BezierCurves_t::Ptr_t BezierCurvesPtr_t;

      BezierCurvesPtr_t project3 (
          const Problem& problem,
          const PathPtr_t& path,
          const value_type& distThr)
      {
        LiegroupSpacePtr_t space (problem.robot()->RnxSOnConfigSpace()->vectorSpacesMerged());
        BezierCurvesPtr_t curves (BezierCurves_t::create (space,
            path->initial(), path->end(),
            path->length()
            ));
           // , path->constraints()));

        const Distance& distance = *problem.distance();
        const PathValidationPtr_t& validation = problem.pathValidation();

        const interval_t& tr = path->timeRange();
        pinocchio::LiegroupElement q1  (space),
                                   q2  (space),
                                   cp1 (space),
                                   cp2 (space),
                                   qt  (space);
        vector_t v1 (path->outputDerivativeSize()),
                 v2 (path->outputDerivativeSize()),
                 vt (path->outputDerivativeSize());
        value_type t1 = tr.first, dt = (tr.second-tr.first)/2;

        q1 = path->initial();
        velocityAt (path, t1, qt.vector(), v1);
        curves->velocityAtControlPoint(0, v1);
        // set velocity at t1
        while (true) {
          bool success = configAndVelocityAt (path, t1+dt, q2.vector(), v2);
          if (!success) {
            hppDout (info, "Could not project configuration at param " << t1+dt);
          }
          value_type dist = distance(q1.vector(), q2.vector());
          cp1 = q1;
          cp1 += (-dist/3)*v1;
          cp2 = q2;
          cp2 += (-dist/3)*v2;

          // Compute hermite distance
          value_type hermite_length = v1.norm() + (cp2 - cp1).norm() + v2.norm();
          if (!success || hermite_length > distThr) {
            hppDout (info, "Distance too big at " << t1 << " (" << dt << ")"
                << " (" << hermite_length << ")");
            dt /= 2;
            continue;
          }

          // Do collision checking
          // TODO Avoid this memory allocation to increase performances.
          BezierCurvesPtr_t p = BezierCurves_t::create (q1, q2, dist);
          p->velocityAtControlPoint(0, v1);
          p->velocityAtControlPoint(1, v2);
          PathPtr_t vp;
          PathValidationReportPtr_t report;
          HPP_START_TIMECOUNTER(validation);
          bool valid = problem.pathValidation ()->validate (p, false, vp, report);
          HPP_STOP_TIMECOUNTER(validation);
          if (!valid) {
            hppDout (info, "q1 = " << setpyformat << one_line(q1.vector())
                << '\n' << "v1 = " << one_line(v1)
                << '\n' << "v2 = " << one_line(v2)
                << '\n' << "q2 = " << one_line(q2.vector())
                << '\n' << "control points\n" << pretty_print(p->controlPoints())
                << unsetpyformat);
            hppDout (info, "Path not valid at " << t1 << " (" << dist << ")\n"
                << *report);
            dt /= 2;
            continue;
          }
          bool first = (t1 == tr.first);
          t1 += dt;
          size_type idx;
          if (t1 >= tr.second)
            idx = (curves->controlPoints().cols()-1)/3;
          else
            idx = curves->insert(t1, q2.vector());
          if (first) { assert (idx==1); }
          curves->velocityAtControlPoint(idx-1, v1);
          curves->velocityAtControlPoint(idx  , v2);
          q1 = q2;
          v1 = v2;
          dt *= 2;
          if (t1 >= tr.second) break;
          if (t1 + dt >= tr.second) dt = tr.second - t1;
        }
        return curves;
      }

      RecursiveHermitePtr_t RecursiveHermite::create (const Problem& problem,
          const value_type& M, const value_type& beta)
      {
        return RecursiveHermitePtr_t (new RecursiveHermite (problem, M, beta));
      }

      RecursiveHermitePtr_t RecursiveHermite::createFromParameters (const Problem& problem)
      {
        value_type beta = problem.getParameter("RecursiveHermite/beta").floatValue();
        if (beta < 0.5 || 1 < beta)
          throw std::invalid_argument ("Parameter \"RecursiveHermite/beta\" "
              "should be between 0.5 and 1");
        value_type M = problem.getParameter("RecursiveHermite/LipschitzConstant").floatValue();
        if (M <= 0)
          throw std::invalid_argument ("Parameter \"RecursiveHermite/LipschitzConstant\" "
              "should greater than 0");
        return create (problem, M, beta);
      }

      RecursiveHermite::RecursiveHermite (const Problem& problem,
          const value_type& M, const value_type& beta) :
        PathOptimizer (problem),
        sm_ (steeringMethod::Hermite::create(problem)),
        M_ (M), beta_ (beta)
      {
        // beta should be between 0.5 and 1.
        if (beta_ < 0.5 || 1 < beta_)
          throw std::invalid_argument ("Beta should be between 0.5 and 1");
      }

      PathVectorPtr_t RecursiveHermite::optimize (const PathVectorPtr_t& input)
      {
        const value_type errThr = problem().getParameter("RecursiveHermite/errorThreshold").floatValue();
        if (errThr <= 0)
          throw std::invalid_argument ("Parameter \"RecursiveHermite/errorThreshold\" "
              "should greater than 0");
        const size_type projectMethod = problem().getParameter("RecursiveHermite/projectionMethod").intValue();

        PathVectorPtr_t flat = PathVector::create
          (input->outputSize(), input->outputDerivativeSize());
        input->flatten(flat);
        PathVectorPtr_t path = PathVector::create
          (input->outputSize(), input->outputDerivativeSize());
        std::vector<PathPtr_t> paths;
        paths.reserve (flat->numberPaths());

        for (std::size_t i = 0; i < flat->numberPaths(); ++i) {
          PathPtr_t p = flat->pathAtRank (i);
          // Remove zero length path.
          if (p->length() > 0)
            paths.push_back(p);
        }
        HermitePaths_t hermites;

        // Make initial vector of hermite curves with continuous velocities
        PathPtr_t prev, next;
        Configuration_t qtmp (input->outputSize());
        vector_t vtmp0 (input->outputDerivativeSize()),
                 vtmp1 (vtmp0.size()),
                 v0 (vtmp0.size()),
                 v1 (vtmp0.size());
        for (std::size_t i = 0; i < paths.size(); ++i)
        {
          PathPtr_t cur = paths[i];
          if (i < paths.size()-1) next = paths[i+1];
          else                    next.reset();

          velocity (prev, cur, qtmp, vtmp0, vtmp1, v0);
          velocity (cur, next, qtmp, vtmp0, vtmp1, v1);

          // TODO handle projected path.
          /*
          InterpolatedPathPtr_t ip = HPP_DYNAMIC_PTR_CAST(InterpolatedPath, path);
          if (ip) {
            typedef InterpolatedPath::InterpolationPoints_t IPs_t;
            const IPs_t& ips = ip->interpolationPoints();
            ps.reserve(ips.size() - 1);
            IPs_t::const_iterator _ip1 = ips.begin(); std::advance (_ip1, 1);
            for (IPs_t::const_iterator _ip0 = ips.begin();
                _ip1 != ips.end(); ++_ip0) {
              ps.push_back (HPP_DYNAMIC_PTR_CAST(Hermite,
                    steer (_ip0->second, _ip1->second)));
              ++_ip1;
            }
            */
          hermites.push_back(sm_->steer (cur->initial(), cur->end(), v0, v1, cur->length()));
          prev = cur;
        }

        // For each segment, do apply the projection algorithm.
        PathVectorPtr_t res = PathVector::create
          (input->outputSize (), input->outputDerivativeSize ());
        for (std::size_t i = 0; i < hermites.size(); ++i) {
          hppDout (info, "Start path " << i);

          value_type thr = infty;
          ConfigProjectorPtr_t proj = getConfigProj (paths[i]);
          // When projector dimension is 0, there are only explicit constraints.
          // TODO At the moment, explicit constraints are not taken into account
          // in the velocity calculation.
          if (proj && proj->dimension() > 0)
            thr = 2 * proj->errorThreshold() / M_;
          if (projectMethod==1)
            res->appendPath(project (problem(), paths[i], thr));
          else if (projectMethod==2)
            res->appendPath(project2(problem(), sm_, paths[i], thr));
          else if (projectMethod==3)
            res->appendPath(project3(problem(), paths[i], thr));
        }

        HPP_DISPLAY_TIMECOUNTER(validation);
	return res;
      }

      bool RecursiveHermite::recurse (const PathPtr_t input,
          const value_type& t0, const value_type& t1,
          const HermitePtr_t& path, PathVectorPtr_t& proj,
          const value_type& acceptThr) const
      {
        if (path->hermiteLength() < acceptThr) {
          // If there are collision, then the path must be split.
          PathPtr_t validPart;
	  PathValidationReportPtr_t report;
          bool valid = problem ().pathValidation ()->validate (path, false, validPart, report);
          if (valid) {
            // TODO this does not work because it is not possible to remove
            // constraints from a path.
            // proj->appendPath (path->copy (ConstraintSetPtr_t()));
            proj->appendPath(path);
            return true;
          }
          return false;
        }

        // Split path into two.
        //const value_type t = 0.5; //path->timeRange().first + path->length() / 2;
        const value_type t = (t0 + t1)/2; //path->timeRange().first + path->length() / 2;
        bool success;
        Configuration_t q1((*input) (t, success));
        if (!success) {
          hppDout (info, "RHP stopped because it could not project a configuration");
          return false;
        }
        const Configuration_t q0 = path->initial ();
        const Configuration_t q2 = path->end ();
        // Velocities must be divided by two because each half is rescale
        // from [0, 0.5] to [0, 1]
        vector_t vHalf (input->outputDerivativeSize());
        velocityAt (input, t, q1, vHalf);

        HermitePtr_t left  = sm_->steer (q0, q1, path->v0(), vHalf, t - t0);
        HermitePtr_t right = sm_->steer (q1, q2, vHalf, path->v1(), t1 - t);

        if (!recurse (input, t0, t, left , proj, acceptThr)) return false;
        if (!recurse (input, t, t1, right, proj, acceptThr)) return false;

        return true;
      }

      // ----------- Declare parameters ------------------------------------- //

      HPP_START_PARAMETER_DECLARATION(pathOptimization_RecursiveHermite)
      Problem::declareParameter(ParameterDescription (Parameter::FLOAT,
            "RecursiveHermite/errorThreshold",
            "The constraints satisfaction threshold used when there are no constraint.",
            Parameter(1.)));
      Problem::declareParameter(ParameterDescription (Parameter::FLOAT,
            "RecursiveHermite/LipschitzConstant",
            "A Lipschitz constant of the constraints.",
            Parameter(10.)));
      Problem::declareParameter(ParameterDescription (Parameter::FLOAT,
            "RecursiveHermite/beta",
            "See \"Fast Interpolation and Time-Optimization on Implicit Contact Submanifolds\" from Kris Hauser.",
            Parameter(0.9)));
      Problem::declareParameter(ParameterDescription (Parameter::INT,
            "RecursiveHermite/projectionMethod",
            "Select the projecction method.",
            Parameter(size_type(1))));
      HPP_END_PARAMETER_DECLARATION(pathOptimization_RecursiveHermite)
    } // namespace pathOptimization
  } // namespace core
} // namespace hpp
