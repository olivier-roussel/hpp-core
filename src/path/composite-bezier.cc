// Copyright (c) 2015, Joseph Mirabel
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

#include <hpp/core/path/composite-bezier.hh>

#include <hpp/util/debug.hh>

#include <hpp/pinocchio/liegroup.hh>
#include <hpp/pinocchio/configuration.hh>

#include <hpp/core/config-projector.hh>
#include <hpp/core/projection-error.hh>
//#include <hpp/core/path/spline.hh>

#include <path/math.hh>

namespace hpp {
  namespace core {
    namespace path {
      using pinocchio::LiegroupElementConstRef;
      using pinocchio::LiegroupElementRef;

      template <int Degree> struct spline_basis_function
      {
        enum { NbCoeffs = Degree + 1 };
        typedef Eigen::Matrix<size_type, NbCoeffs, 1> Factorials_t;
        typedef Eigen::Matrix<value_type, NbCoeffs, 1> Coeffs_t;
        typedef Eigen::Matrix<value_type, NbCoeffs, NbCoeffs> IntegralCoeffs_t;

        static void derivative (const size_type order, const value_type& t, Coeffs_t& res);
        static void velocity_bound (const value_type& t0, const value_type& t1, Coeffs_t& res);

        private:
        static Coeffs_t absBound (bool up);
      };
      template <int Degree>
      void spline_basis_function<Degree>::derivative
      (const size_type k, const value_type& t, Coeffs_t& res)
      {
        res.setZero();
        if (k > Degree) return;
        static Factorials_t factors = binomials<NbCoeffs>::factorials();
        Coeffs_t powersOfT, powersOfOneMinusT;
        powersOfT        (0) = 1;
        powersOfOneMinusT(0) = 1;
        const value_type oneMinusT = 1 - t;
        for (size_type i = 1; i < NbCoeffs; ++i) {
          powersOfT        (i) = powersOfT        (i - 1) * t;
          powersOfOneMinusT(i) = powersOfOneMinusT(i - 1) * oneMinusT;
        }

        for (size_type i = 0; i < NbCoeffs; ++i) {
          for (size_type p = std::max((size_type)0, k + i - Degree); p <= std::min(i, k); ++p) {
            size_type ip = i - p;
            res(i) += value_type (( (k - p) % 2 == 0 ? 1 : -1 )
              * binomials<NbCoeffs>::binomial (k, p))
              *   powersOfT(ip) * powersOfOneMinusT(Degree - k - ip)
              / value_type( factors  (ip) * factors  (Degree - k - ip) );
          }
        }
        res *= value_type(factors(Degree));
      }
      template <int Degree>
      void spline_basis_function<Degree>::velocity_bound
      (const value_type& t0, const value_type& t1, Coeffs_t& res)
      {
        static const Coeffs_t b_up = absBound(true);
        static const Coeffs_t b_um = absBound(false);
        Coeffs_t bt0, bt1;
        derivative(1, t0, bt0);
        derivative(1, t1, bt1);

        res.noalias() = bt0.cwiseAbs().cwiseMax(bt1.cwiseAbs());

        // Case i = 0 and i = n
        // Nothing to do.

        // Case i = 1, n-1
        if (t0 * Degree < 2 && 2 < t1 * Degree) // If max for i = 1 in [t0, t1]
          res(1) = std::max(res(1), b_up(1));
        if (t0 * Degree < Degree - 2 && Degree - 2 < t1 * Degree) // If max for i = n-1 in [t0, t1]
          res(Degree-1) = std::max(res(Degree-1), b_up(Degree-1));

        // Case 2 <= i <= n-2
        if (Degree > 3) {
          // when u_p(i) in [t0, t1], consider b_up.
          size_type r1 = std::max(size_type(std::ceil (t0 * (Degree - 1))), size_type(2)),
                    r2 =          size_type(std::floor(t1 * (Degree - 1))),
                    nr = std::max(r2 - r1 + 1, size_type(0));
          res.segment(r1, nr).noalias() = res.segment(r1, nr).cwiseMax(b_up.segment(r1, nr));

          // when u_m(i) in [t0, t1], consider b_um.
          r1 =          size_type(std::ceil (1 + t0 * (Degree - 1))),
          r2 = std::min(size_type(std::floor(1 + t1 * (Degree - 1))), size_type(Degree - 2));
          nr = std::max(r2 - r1 + 1, size_type(0));
          res.segment(r1, nr).noalias() = res.segment(r1, nr).cwiseMax(b_um.segment(r1, nr));
        }
      }
      template <int Degree>
      typename spline_basis_function<Degree>::Coeffs_t spline_basis_function<Degree>::absBound (bool up)
      {
        Coeffs_t res;
        Factorials_t iToThePowerOfI (Factorials_t::Ones());
        for (size_type i = 1; i < Degree; ++i) {
          for (size_type j = 0; j < i; ++j)
            iToThePowerOfI(i) *= i;
        }
        res(0) = res(Degree) = Degree;
        if (Degree > 1) {
          if (Degree == 2) res(1) = res(Degree - 1) = 2;
          else res(1) = res(Degree - 1)
            = value_type(iToThePowerOfI(Degree - 2)) / std::pow (Degree, Degree-3);
          for (size_type i = 2; i < Degree - 1; ++i) {
            const size_type p = (up
                ? iToThePowerOfI(i) * iToThePowerOfI(Degree - i - 1)
                : iToThePowerOfI(Degree - i) * iToThePowerOfI(i - 1));
            res(i) = value_type(binomials<NbCoeffs>::binomial(Degree, i) * p) / value_type(iToThePowerOfI(Degree - 1));
          }
        }
        return res;
      }

      template <int Order>
      CompositeBezier<Order>::CompositeBezier (
          const LiegroupSpacePtr_t& space,
          ConfigurationIn_t init,
          ConfigurationIn_t end,
          value_type length) :
        Path (interval_t (0, length), space->nq (), space->nv ()),
        space_ (space),
        params_(2),
        controlPoints_(space->nq(), Order+1),
        npts_(Order+1)
      {
        assert (space);
        assert (length >= 0);
        assert (init.size() == space->nq ());
        assert (end .size() == space->nq ());

        params_[0] = 0;
        params_[1] = length;

        controlPoints_. leftCols<(Order+1)/2>() = init.replicate<1,(Order+1)/2>();
        controlPoints_.rightCols<(Order+1)/2>() =  end.replicate<1,(Order+1)/2>();

        assert (!constraints ());
      }

      template <int Order>
      CompositeBezier<Order>::CompositeBezier (const LiegroupSpacePtr_t& space,
                                  ConfigurationIn_t init,
                                  ConfigurationIn_t end,
                                  value_type length,
                                  ConstraintSetPtr_t constraints) :
        Path (interval_t (0, length), space->nq (), space->nv (), constraints),
        space_ (space),
        params_(2),
        controlPoints_(space->nq(), Order+1),
        npts_(Order+1)
      {
        assert (space);
        assert (length >= 0);
        assert (init.size() == space->nq ());
        assert (end .size() == space->nq ());

        params_[0] = 0;
        params_[1] = length;

        controlPoints_. leftCols<(Order+1)/2>() = init.replicate<1,(Order+1)/2>();
        controlPoints_.rightCols<(Order+1)/2>() =  end.replicate<1,(Order+1)/2>();
      }

      template <int Order>
      CompositeBezier<Order>::CompositeBezier (const CompositeBezier& path) :
        Path (path),
        space_ (path.space_),
        params_ (path.params_),
        controlPoints_ (path.controlPoints_),
        npts_ (path.npts_)
      {}

      template <int Order>
      CompositeBezier<Order>::CompositeBezier (const CompositeBezier& path,
                                  const ConstraintSetPtr_t& constraints) :
        Path (path, constraints),
        space_ (path.space_),
        params_ (path.params_),
        controlPoints_ (path.controlPoints_),
        npts_ (path.npts_)
      {}

      template <int Order>
      void CompositeBezier<Order>::init (Ptr_t self)
      {
        Path::init (self);
        weak_ = self;
        checkPath ();
      }

      template <int Order>
      void CompositeBezier<Order>::initCopy (Ptr_t self)
      {
        Path::init (self);
        weak_ = self;
        checkPath ();
      }

      /// \param t time in [0,1].
      template <typename ControlPoints>
      void recursiveBezier (
          const LiegroupSpacePtr_t& space,
          const value_type& t,
          const Eigen::MatrixBase<ControlPoints>& pts, size_type start, size_type N,
          matrix_t& qs, size_type i)
      {
        assert (0 <= t && t <= 1);

        matrix_t::ColXpr q = qs.col(i);
        if (N == 1) {
          // This should not happen
          // Return pts.col(0)
          q.noalias() = pts.derived().col(start);
          return;
        } else if (N == 2) {
          LiegroupElementConstRef A0 (pts.derived().col(start  ), space),
                                  A1 (pts.derived().col(start+1), space);
          // qs = A0 + t * (A1-A0)
          q = (A0 + t * (A1-A0)).vector();
        } else {
          recursiveBezier(space, t,
              pts, start+1, N-1,
              qs, i  );
          recursiveBezier(space, t,
              pts, start  , N-1,
              qs, i+1);
          LiegroupElementConstRef A0 (qs.col(i+1), space),
                                  A1 (qs.col(i  ), space);
          // qs = A0 + t * (A1-A0)
          q = (A0 + t * (A1-A0)).vector();
        }
      }

      /// \param t time in [0,1].
      template <typename ControlPoints>
      void recursiveJBezier (
          const LiegroupSpacePtr_t& space,
          const value_type& t,
          const Eigen::MatrixBase<ControlPoints>& pts, size_type start, size_type N,
          matrix_t& qs, size_type i,
          matrix_t& vs, size_type j)
      {
        //std::cout << start << " " << N << " " << t << ":" << incindent << iendl;
        using pinocchio::DerivativeTimesInput;

        assert (0 <= t && t <= 1);

        matrix_t::ColXpr q = qs.col(i);
        matrix_t::ColXpr v = vs.col(j);
        if (N == 1) {
          // This should not happen
          q = pts.derived().col(start);
          v.setZero();
        } else if (N == 2) {
          // res = d+/dv * (A1-A0)
          LiegroupElementConstRef A0 (pts.derived().col(start  ), space),
                                  A1 (pts.derived().col(start+1), space);
          matrix_t::ColXpr w = vs.col(j+1);
          v = A1-A0;
          w = t*v;
          space->dIntegrate_dv<DerivativeTimesInput> (A0, w, v);
          q = (A0 + w).vector();
          //std::cout << pts.template middleCols<2>(start) << iendl;
        } else {
          // qs.col(i+1) <- A0
          // vs.col(j+1) <- JA0
          recursiveJBezier(space, t,
              pts, start, N-1,
              qs, i+1,
              vs, j+1);
          // qs.col(i+2) <- A1
          // vs.col(j+2) <- JA1
          recursiveJBezier(space, t,
              pts, start+1, N-1,
              qs, i+2,
              vs, j+2);

          LiegroupElementConstRef A0 (qs.col(i+1), space),
                                  A1 (qs.col(i+2), space);
          matrix_t::ColXpr a0 (qs.col(i+1));
          matrix_t::ColXpr a1 (qs.col(i+2));
          matrix_t::ColXpr JA0   (vs.col(j+1));
          matrix_t::ColXpr JA1   (vs.col(j+2));
          matrix_t::ColXpr w     (vs.col(j+3));
          matrix_t::ColXpr dv_dt (vs.col(j+4));

          // res = A0 + t * (A1-A0)
          // TODO might not be always needed.
          v = A1-A0;
          w = t * v;
          q = (A0 + w).vector();
          //std::cout << "v " << v << iendl;

          // dv/dt = d-/dq1 * JA1 + d-/dq0 * JA0
          dv_dt = JA0;
          space->dDifference_dq0<DerivativeTimesInput> (a0, a1, dv_dt);
          space->dDifference_dq1<DerivativeTimesInput> (a0, a1, JA1);
          dv_dt.noalias() += JA1;
          //std::cout << "dv_dt " << dv_dt << iendl;

          // d+/dq JA0 + d+/dv * (t * dv/dt + (A1-A0))
          v += t*dv_dt;
          space->dIntegrate_dv<DerivativeTimesInput> (A0, w, v);
          space->dIntegrate_dq<DerivativeTimesInput> (A0, w, JA0);
          v += JA0;

        }
        //std::cout
          //<< q.transpose() << iendl
          //<< v.transpose() << decindent << iendl;
      }

      template <int Order>
      bool CompositeBezier<Order>::impl_compute (ConfigurationOut_t result,
                                       value_type param) const
      {
        assert (param >= paramRange().first);
        if (param == paramRange ().first || paramLength() == 0) {
          result.noalias () = controlPoints_.col(0);
          return true;
        }
        if (param >= params_.back()) {
          result.noalias () = controlPoints_.col(npts_-1);
          return true;
        }

        value_type localParam, T;
        const size_type idx = indexAtParam (param, T, localParam);
        const value_type u = localParam / T;

        matrix_t qs (space_->nq(), Order+2);
        recursiveBezier (space_, u,
            controlPoints_.middleCols<Order+1>(idx), 0, Order+1,
            qs, 0);
        result = qs.col(0);
        return true;
      }

      template <int Order>
      void CompositeBezier<Order>::impl_derivative
      (vectorOut_t result, const value_type& param, size_type order) const
      {
        if (   order > Order
            || paramRange ().first == paramRange ().second
            ) {
          result.setZero ();
          return;
        }

        value_type localParam, T;
        const size_type idx = indexAtParam (param, T, localParam);
        const value_type u = localParam / T;

        assert (param >= paramRange().first);

        if (order == 1) {
          matrix_t qs (space_->nq(), Order+3);
          matrix_t vs (space_->nv(), Order+5);

          recursiveJBezier (space_, u,
              controlPoints_.middleCols<Order+1>(idx), 0, Order+1,
              qs, 0, vs, 0);
          result.noalias() = vs.col(0) / T;
        }
      }

      template <int Order>
      void CompositeBezier<Order>::impl_velocityBound (vectorOut_t result,
          const value_type& param0, const value_type& param1) const
      {
        throw std::logic_error ("unimplemented");
        value_type localParam0, T0,
                   localParam1, T1;
        const size_type idx0 = indexAtParam (param0, T0, localParam0),
                        idx1 = indexAtParam (param1, T1, localParam1);

        result.setZero();

        typedef spline_basis_function<Order+1> sbf_t;
        //Eigen::Matrix<value_type, Eigen::Dynamic, Order+1> vs (space_->nv(), Order+1);
        //vs.col(0).setZero();
        Eigen::Matrix<value_type, Eigen::Dynamic, Order> vs (space_->nv(), Order);

        typename sbf_t::Coeffs_t ub;
        value_type u0 = localParam0/T0, u1 = localParam1/T1;
        for (size_type i = idx0; i < idx1+1; i+=Order) {
          // Velocities
          LiegroupElementConstRef A0 (controlPoints_.col(i), space_);
          for (size_type k = 1; k < vs.cols(); ++k)
            vs.col(k-1) = (LiegroupElementConstRef(controlPoints_.col(i+k), space_) - A0).cwiseAbs();
          // Velocity bound
          if (i==idx0) {
            sbf_t::velocity_bound (u0, (idx0==idx1?u1:1), ub);
            u0 = 0;
            ub /= T0;
          } else if (i==idx1) {
            sbf_t::velocity_bound (0, u1, ub);
            ub /= T1;
          } else {
            value_type T = params_[i] - params_[i-1];
            ub.setConstant (2*(Order+1)/T);
          }
          result.noalias() = result.cwiseMax(vs * ub.template tail<Order>());
        }
      }

      template <int Order>
      template <typename Derived>
      CompositeBezier<Order>::CompositeBezier (const LiegroupSpacePtr_t& space,
            const Params_t& params,
            const Eigen::DenseBase<Derived>& controlPoints,
            ConstraintSetPtr_t constraints) :
        Path (interval_t (0, params_.back()), space->nq (), space->nv (), constraints),
        space_ (space),
        params_ (params),
        controlPoints_ (controlPoints.derived()),
        npts_ (controlPoints.derived().cols())
      {}

      template <int Order>
      PathPtr_t CompositeBezier<Order>::reverse () const
      {
        Params_t params (params_.size());
        params[0] = 0;
        for (std::size_t i = 1; i < params.size(); ++i)
          params[i] = params[i-1] +
            params_[params_.size()-i] - params_[params_.size()-i-1];
        Ptr_t result (new CompositeBezier(space_, params,
              controlPoints_.leftCols(npts_).rowwise().reverse(),
              constraints()));
        result->init (result);
        return result;
      }

      template <int Order>
      size_type CompositeBezier<Order>::insert (const value_type& param, ConfigurationIn_t config)
      {
        assert((npts_-1)%Order == 0);
        value_type T, localParam;
        const size_type idx = indexAtParam (param, T, localParam);
        assert (idx >= 0 && idx < params_.size());
        params_.insert (params_.begin()+idx+1, param);
        // We want to insert Order control points
        // after index idx+(Order+1)/2.
        const size_type start = idx + (Order+1)/2;

        // Do we need reallocation
        if (controlPoints_.cols() < npts_+Order) {
          ControlPoints_t ctrlpts (space_->nq(), (npts_+Order)*2);
          // Copy from 0 to start.
          ctrlpts.leftCols(start) = controlPoints_.leftCols(start);
          // Copy from start to npts_, order column shifting to the right.
          ctrlpts           .middleCols(start+Order, npts_-start)
            = controlPoints_.middleCols(start      , npts_-start);
          controlPoints_.swap(ctrlpts);
        } else {
          // Order column shift to the right of the columns from start.
          for (size_type c = npts_-1; c > start; c -= Order)
            controlPoints_    .middleCols(c+Order, Order).noalias()
              = controlPoints_.middleCols(c      , Order);
        }

        controlPoints_.middleCols<Order>(start) = config.replicate<1,Order>();
        npts_ += Order;
        return (idx/Order)+1;
      }

      template <int Order>
      void CompositeBezier<Order>::velocityAtControlPoint (const size_type& iPt, vectorIn_t v)
      {
        assert(Order>=3);
        assert (iPt*Order < npts_);

        const size_type idx = iPt*Order;
        LiegroupElementConstRef q (controlPoints_.col(idx  ), space_);
        vector_t dq (v.size()); 
        if (iPt > 0) {
          LiegroupElementRef qprev (controlPoints_.col(idx-1), space_);
          const value_type T = params_[iPt] - params_[iPt-1];
          dq.noalias() = (-T/3)*v;
          qprev = q;
          qprev += dq;
        }
        if (idx < npts_-1) {
          LiegroupElementRef qnext (controlPoints_.col(idx+1), space_);
          const value_type T = params_[iPt+1] - params_[iPt];
          dq.noalias() = (T/3)*v;
          qnext = q;
          qnext += dq;
        }
      }

      template <int Order>
      size_type CompositeBezier<Order>::indexAtParam (const value_type& param,
          value_type& T,
          value_type& localParam) const
      {
        assert(param >= 0 && param <= paramLength());

        // *(_param-1) <= param < *_param
        Params_t::const_iterator _param =
          std::upper_bound (params_.begin(), params_.end(), param);

        std::size_t rank = _param - params_.begin();
        if (_param == params_.end()) {
          rank = params_.size()-1;
        }

        assert (rank>0);
        T =          params_[rank] - params_[rank-1];
        localParam = param         - params_[rank-1];
        assert(localParam<=T);
        return (rank-1) * Order;
      }

      template class CompositeBezier<1>;
      template class CompositeBezier<3>;
    } //   namespace path
  } //   namespace core
} // namespace hpp

