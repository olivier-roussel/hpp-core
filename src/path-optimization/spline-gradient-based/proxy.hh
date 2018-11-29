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

#ifndef SRC_PATH_OPTIMIZATION_SPLINE_GRADIENT_BASED_PROXY_HH
#define SRC_PATH_OPTIMIZATION_SPLINE_GRADIENT_BASED_PROXY_HH

#include <hpp/constraints/differentiable-function.hh>

#include <hpp/core/config.hh>

namespace hpp {
  namespace core {
    namespace pathOptimization {
      /// This class computes an output spline from an input spline.
      /// Let \f$ I = P^a_i \f$ be the input spline of lenth T_i
      /// and \f$ t \in [0,T_i]\f$.
      ///
      /// This class computes an output function \f$ O = P^b_i \f$, whose input
      /// variables are the parameters of spline \f$ I \f$, such that:
      /// - \f$ O(t) = inner (I(t)) \f$
      /// - \f$ \frac{dO}{dt}(t) = Jinner (I(t)) \frac{dI}{dt}(t) \f$
      template <int _PolynomeBasis, int _SplineOrder>
      class HPP_CORE_LOCAL ProxyFunction :
        public constraints::DifferentiableFunction
      {
        public:
          typedef boost::shared_ptr<ProxyFunction> Ptr_t;
          typedef path::Spline<_PolynomeBasis, _SplineOrder> Spline;

          /// Constructor with no right hand side
          /// \param base spline base
          /// \param u ratio at which this function applies.
          /// \param T length of the spline. (May be one day this will be an optimization variable.)
          ProxyFunction (const DifferentiableFunctionPtr_t& inner,
              const LiegroupElement& ba_i,
              const value_type& u,
              const value_type& T_i,
              const pinocchio::ArrayXb& rhs) :
            DifferentiableFunction (
                (rhs.count() + inner->inputDerivativeSize())*Spline::NbCoeffs,
                (rhs.count() + inner->inputDerivativeSize())*Spline::NbCoeffs,
                inner->outputSpace(), "Proxy on " + inner->name()),
            inner_ (inner),
            ba_i_ (ba_i),
            u_ (u),
            T_i_ (T_i),
            rhs_ (rhs),
            Nr_ (rhs.count()                 *Spline::NbCoeffs),
            Ni_ (inner->inputDerivativeSize()*Spline::NbCoeffs),
            bfv_ (computeBFV(u_))
          {
            activeParameters_          .setConstant(true);
            activeDerivativeParameters_.setConstant(true);
          }

          /// Constructor with a right hand side
          ProxyFunction (const DifferentiableFunctionPtr_t& inner,
              const LiegroupElement& ba_i,
              const value_type& u,
              const value_type& T_i,
              const DifferentiableFunctionPtr_t& rhs,
              const segment_t& innerIn,
              const segment_t& rhsIn) :
            DifferentiableFunction (
                (rhs.count() + inner->inputDerivativeSize())*Spline::NbCoeffs,
                (rhs.count() + inner->inputDerivativeSize())*Spline::NbCoeffs,
                inner->outputSpace(), "Proxy on " + inner->name()),
            inner_ (inner),
            ba_i_ (ba_i),
            u_ (u),
            T_i_ (T_i),
            rhs_ (rhs),
            Nr_ (rhs.count()                 *Spline::NbCoeffs),
            Ni_ (inner->inputDerivativeSize()*Spline::NbCoeffs),
            bfv_ (computeBFV(u_))
          {
            activeParameters_          .setConstant(true);
            activeDerivativeParameters_.setConstant(true);
          }

        protected:
          typedef Eigen::Map<Spline::ParameterMatrix_t> ParameterMap_t;

          void spline (const vectorIn_t parameters, Configuration_t& q) const
          {
            Eigen::Map<const Spline::ParameterMatrix_t> P_i (parameters);
            // TODO Temporary check to make sure storage order is correct.
            assert (P_i (0, 1) == parameters(1));

            q.resize   (ba_i_.nq());
            vector_t v (ba_i_.nv());
            Spline::value (ba_i_, P_i, u_, q, v);
          }

          void impl_compute (LiegroupElement& y, vectorIn_t parameters) const
          {
            Configuration_t q;
            spline (parameters.segment(inner_.first,inner_.second), q);
            inner_->value(y, q);

            if (rhs_) {
            }
          }

          void impl_jacobian (matrixOut_t J, vectorIn_t parameters) const
          {
            Configuration_t q;
            spline (parameters, q);

            // matrix_t Jinner (outputDerivativeSize(), inner_->inputDerivativeSize());
            const size_type Nv = inner_->inputDerivativeSize();

            // Store Jinner in the first Nv columns
            matrixOut_t Jinner (J.leftCols(Nv));
            inner_->jacobian(Jinner, q);

            // Compute df/dPi = paramDerCoeffs(i) * Jinner
            // Because Jinner is store at the location of df/dP0
            // we must compute the value to i > 0 first.
            size_type i_Pi = Nv;
            for (size_type i = 1; i < Spline::NbCoeffs; ++i) {
              J.middleCols (i_Pi, Nv).noalias() = bfv(i) * Jinner;
              i_Pi += Nv;
            }
            assert (inputDerivativeSize() == i_Pi);
            // Compute df/dP0
            Jinner *= bfv(0);
          }

          std::ostream& print (std::ostream& os) const
          {
            constraints::DifferentiableFunction::print(os);
            os << "at ratio " << u_;
          }

          static Spline::BasisFunctionVector_t computeBFV (const value_type u)
          {
            Spline::BasisFunctionVector_t bfv;
            Spline::timeFreeBasisFunctionDerivative (0, u, bfv);
            return bfv;
          }

          DifferentiableFunctionPtr_t inner_;
          LiegroupElement ba_i_;
          const value_type u_, T_i_;

          DifferentiableFunctionPtr_t rhs_;
          const segment_t innerIn_, rhsIn_;

          const Spline::BasisFunctionVector_t bfv_;
      }; // class ProxyFunction

      class ProxiedFunction :
        public constraints::DifferentiableFunction
      {
      };
    } // namespace pathOptimization
  } // namespace core
} // namespace hpp

#endif // SRC_PATH_OPTIMIZATION_SPLINE_GRADIENT_BASED_PROXY_HH
