// Copyright (c) 2018, Joseph Mirabel
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

#include <hpp/core/path-optimization/quadratic-program.hh>

#include <hpp/util/timer.hh>

#include <path-optimization/spline-gradient-based/eiquadprog-fast.hpp>

namespace hpp {
  namespace core {
    /// \addtogroup path_optimization
    /// \{
    namespace pathOptimization {
      HPP_DEFINE_TIMECOUNTER(QuadraticProgram_decompose);
      HPP_DEFINE_TIMECOUNTER(QuadraticProgram_computeLLT);
      HPP_DEFINE_TIMECOUNTER(QuadraticProgram_solve_quadprog);

      void QuadraticProgram::initQuadProg ()
      {
        quadprog_ = new EiquadprogFast ();
      }

      QuadraticProgram::~QuadraticProgram ()
      {
        HPP_DISPLAY_TIMECOUNTER(QuadraticProgram_decompose);
        HPP_DISPLAY_TIMECOUNTER(QuadraticProgram_computeLLT);
        HPP_DISPLAY_TIMECOUNTER(QuadraticProgram_solve_quadprog);
        if (quadprog_) delete quadprog_;
      }

      void QuadraticProgram::decompose ()
      {
        HPP_SCOPE_TIMECOUNTER(QuadraticProgram_decompose);
        dec.compute(H);
        assert(dec.rank() == H.rows());
      }

      void QuadraticProgram::computeLLT()
      {
        HPP_SCOPE_TIMECOUNTER(QuadraticProgram_computeLLT);
        trace = H.trace();
        llt.compute(H);
      }

      double QuadraticProgram::solve(const LinearConstraint& ce, const LinearConstraint& ci)
      {
        HPP_SCOPE_TIMECOUNTER(QuadraticProgram_solve_quadprog);

        quadprog_->m_J.setIdentity(H.rows(), H.rows());
        llt.matrixU().solveInPlace(quadprog_->m_J);
        quadprog_->is_inverse_provided_ = true;

        // min   0.5 * x G x + g0 x
        // s.t.  CE^T x + ce0 = 0
        //       CI^T x + ci0 >= 0
        EiquadprogFast_status status = quadprog_->solve_quadprog (
            H, b,
            ce.J.transpose(), - ce.b,
            ci.J.transpose(), - ci.b,
            xStar);
        activeSetSize = (int)quadprog_->getActiveSetSize ();
        activeConstraint = quadprog_->getActiveSet ();
        switch (status) {
          case EIQUADPROG_FAST_UNBOUNDED:
            hppDout (warning, "Quadratic problem is not bounded");
          case EIQUADPROG_FAST_OPTIMAL:
            break;
          case EIQUADPROG_FAST_INFEASIBLE:
            hppDout (error, "Quadratic problem is not feasible");
            break;
          case EIQUADPROG_FAST_MAX_ITER_REACHED:
            hppDout (error, "Quadratic problem resolution reached the maximum number of iterations.");
            break;
          case EIQUADPROG_FAST_REDUNDANT_EQUALITIES:
            hppDout (error, "Constraint of quadratic problem are linearly dependent.");
            break;
        }
        return quadprog_->getObjValue();
      }
    } // namespace pathOptimization
  }  // namespace core
} // namespace hpp
