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

#ifndef SRC_PATH_OPTIMIZATION_SPLINE_GRADIENT_BASED_CONTINUITY_HH
#define SRC_PATH_OPTIMIZATION_SPLINE_GRADIENT_BASED_CONTINUITY_HH

#include <hpp/constraints/differentiable-function.hh>

#include <hpp/core/config.hh>

namespace hpp {
  namespace core {
    namespace pathOptimization {
      /**
       * Explicit continuity constraint.
       * 
       *  For \f$ 0 \le i < N \f$:
       *  - \f$      P (b_i, c_i, T_i, T_i) =      P (b_{i+1}, c_{i+1}, T_{i+1}, 0) \f$
       *  - \f$ \dot{P}(b_i, c_i, T_i, T_i) = \dot{P}(b_{i+1}, c_{i+1}, T_{i+1}, 0) \f$
       *  - ...
       *  depending on the value of \c order.
       * 
       *  If \f$ i = -1 \f$, then:
       *  - \f$ base_L =      P (b_0, c_0, T_0, 0) \f$
       *  - \f$      0 = \dot{P}(b_0, c_0, T_0, 0) \f$
       *  - ...
       * 
       *  If \f$ i =  N \f$, then:
       *  - \f$      P (b_N, c_N, T_N, T_N) = base_R \f$
       *  - \f$ \dot{P}(b_N, c_N, T_N, T_N) = 0      \f$
       *  - ...
       *
       * \todo T_i and T_{i+1} could be parameters of the function (for non-linear optimization).
       * \todo Add implicit constraint ? It will likely be very slow.
      **/
      template <int _PolynomeBasis, int _SplineOrder>
      class HPP_CORE_LOCAL ExplicitContinuityFunction;

      /// Implementation for Bezier curves.
      template <int _SplineOrder>
      class HPP_CORE_LOCAL ExplicitContinuityFunction <path::BernsteinBasis, _SplineOrder>
        public constraints::DifferentiableFunction
      {
        private:
          enum State { FIRST, MIDDLE, LAST };
        public:
          enum { MaxContinuityOrder = int( (SplineOrder - 1) / 2) };

          typedef boost::shared_ptr<ContinuityFunction> Ptr_t;
          typedef path::Spline<path::BernsteinBasis, _SplineOrder> Spline;

          ContinuityFunction (const DevicePtr_t robot,
              const LiegroupElement& b_i,
              const LiegroupElement& b_i_plus_1,
              const value_type& T_i,
              const value_type& T_i_plus_1,
              const size_type i,
              const size_type N,
              const size_type order = MaxContinuityOrder) :
            DifferentiableFunction (
                ( i == -1 || i == N ) ? 0 : robot->numberDof() * (order + 1), // nb input parameters (arg)
                ( i == -1 || i == N ) ? 0 : robot->numberDof() * (order + 1), // nb input parameters (der)
                LiegroupSpace::create(robot->numberDof() * (order + 1)),
                (std::ostringstream() << "Continuity " << i).str()),
            state_ ( (i == -1) ? FIRST : ( (i==N) ? LAST : MIDDLE)),
            order_ (order),
            b_i_   (b_i),
            b_ip1_ (b_i_plus_1),
            T_i_   (T_i)
            T_ip1_ (T_i_plus_1)
          {
            assert (0 <= order && order <= MaxContinuityOrder);

            activeParameters_          .setConstant(true);
            activeDerivativeParameters_.setConstant(true);

            if (state_ == FIRST || state_ == LAST) v_ = Vector_t::Zero(robot->numberDof());
          }

        protected:
          void impl_compute (LiegroupElement& y, vectorIn_t v) const
          {
            const size_type nv = b_i_.space()->nv();

            // Computes v^1_{i+1} or v^n_i
            switch (state_) {
              case FIRST: // v^1_{i+1}
                y.vector().head (nv) = (b_i_ - b_ip1_).vector();
                break;
              case MIDDLE: // v^n_i
                y.vector().tail (nv) = ((b_ip1_ + v.head(nv)) - b_i_).vector();
                break;
              case LAST: // v^n_i
                y.vector().tail (nv) = (b_ip1_ - b_i_).vector();
                break;
            }

            if (order == 0) return;
            const size_type M = y.size();
            size_type idx;

            // Computes v^2_{i+1} or v^{n-1}_i
            switch (state_) {
              case FIRST: // v^2_{i+1}
                idx = nv;
                y.vector().segment (idx, nv) = (Spline::NbCoeffs / T_ip1_) * v_;
                y.JDifference <true>(
                    b_ip1_.vector(),
                    (b_ip1_ + y.vector().head(nv)).vector(),
                    matrix_t(),
                    y.vector().segment (idx, nv));
                y.vector().segment (idx, nv) += y.vector().head(nv);
                break;
              case MIDDLE: // v^{n-1}_i
                idx = M - 2*nv;
                y.vector().segment(idx, nv).noalias() =
                  (- T_i_ / T_ip1_ ) * (v.segment(nv,nv) - v.head(nv));
                // y <- dI_dq * y
                y.dIntegrate_dq (b_ip1_, v.head(nv), 
                    y.vector().segment (idx, nv));
                // y <- dI_dq^-1 * y
                y.JDifference <true>(
                    b_i_.vector(),
                    (b_i_ + y.vector().tail(nv)).vector(),
                    matrix_t(),
                    y.vector().segment (idx, nv));
                y.vector().segment (idx, nv) += y.vector().tail(nv);
                break;
              case LAST: // v^{n-1}_i
                idx = M - 2*nv;
                y.vector().segment(idx, nv).noalias() = 
                  (- T_i_ / Spline::NbCoeffs) * v_;
                // y <- dI_dq^-1 * y
                y.JDifference <true>(
                    b_i_.vector(),
                    (b_i_ + y.vector().tail(nv)).vector(),
                    matrix_t(),
                    y.vector().segment (idx, nv));
                y.vector().segment (idx, nv) += y.vector().tail(nv);
                break;
            }

            assert (order == 1 && "Currently, higher order of continuity is not supported "
                "because we need d2I/dv2.");
          }

          void impl_jacobian (matrixOut_t J, vectorIn_t parameters) const
          {
            // TODO this works:
            // - for vector spaces,
            // - for SO(3) if b_i_ and b_ip1_ are equal.
            J.setIdentity();
          }

          std::ostream& print (std::ostream& os) const
          {
            constraints::DifferentiableFunction::print(os);
            os << "order " << order_ << ", state " << state_ << iendl;
          }

        private:
          const State state_;
          const size_type order_;
          const LiegroupElement b_i_, b_ip1_;
          const value_type T_i_, T_ip1_;
          vector_t v_;
      }; // class ContinuityFunction
    } // namespace pathOptimization
  } // namespace core
} // namespace hpp

#endif // SRC_PATH_OPTIMIZATION_SPLINE_GRADIENT_BASED_CONTINUITY_HH
