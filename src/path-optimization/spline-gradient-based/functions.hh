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

#include <hpp/constraints/differentiable-function.hh>

#include <hpp/core/config.hh>

namespace hpp {
  namespace core {
    namespace pathOptimization {
      typedef core::segment_t segment_t;
      namespace {
        std::string toStr (const segment_t& s)
        {
          std::ostringstream os;
          os << "[ " << s.first << ", " << s.first + s.second << " ]";
          return os.str();
        }
      }

      /// Apply the constraint on a subspace of the input space.
      /// i.e.: \f$ f (q_0, ... , q_n) = f_{inner} (q_k) \f$
      template <int _PolynomeBasis, int _SplineOrder>
      class HPP_CORE_LOCAL SplineFunction :
        public constraints::DifferentiableFunction
      {
        public:
          typedef boost::shared_ptr<SplineFunction> Ptr_t;
          typedef path::Spline<_PolynomeBasis, _SplineOrder> Spline;

          /// \param base spline base
          /// \param u ratio at which this function applies.
          /// \param T length of the spline. (May be one day this will be an optimization variable.)
          SplineFunction (const DifferentiableFunctionPtr_t& inner,
              const DevicePtr_t robot,
              const Configuration_t& base,
              const value_type& u,
              const value_type& T) :
            DifferentiableFunction (
                inner->inputDerivativeSize()*Spline::NbCoeffs,
                inner->inputDerivativeSize()*Spline::NbCoeffs,
                inner->outputSpace(), inner->name() + " | " + toStr(inArgs)),
            inner_ (inner),
            // base_ (base, robot->configSpace()->vectorSpaceMerged())
            u_ (u),
            spline_ (Spline::create (robot, interval_t(0, T), ConstraintSetPtr_t()))
          {
            // TODO Compute activeParameters_ and activeDerivativeParameters_
            activeParameters_          .setConstant(true);
            activeDerivativeParameters_.setConstant(true);

            // Spline::basisFunctionDerivative(0, u, basisFunc);
            spline_->base(base_);
          }

        protected:
          typedef Eigen::Map<Spline::ParameterMatrix_t> ParameterMap_t;

          void impl_compute (LiegroupElement& y, vectorIn_t parameters) const
          {
            spline_->rowParameters (parameters);
            Configuration_t q (spline_->outputSize());
            (*spline_) (q, u*spline_->length());
            inner_->value(y, q);
          }

          void impl_jacobian (matrixOut_t J, vectorIn_t parameters) const
          {
            spline_->rowParameters (parameters);
            Configuration_t q (spline_->outputSize());
            const value_type t = u*spline_->length();
            (*spline_) (q, t);

            Eigen::Matrix<value_type, Spline::NbCoeffs, 1> paramDerCoeffs;
            spline_->parameterDerivativeCoefficients (paramDerCoeffs, t);

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
              J.middleCols (i_Pi, Nv).noalias() = paramDerCoeffs(i) * Jinner;
              i_Pi += Nv;
            }
            assert (inputDerivativeSize() == i_Pi);
            // Compute df/dP0
            Jinner *= paramDerCoeffs(0);
          }

          std::ostream& print (std::ostream& os) const
          {
            constraints::DifferentiableFunction::print(os);
            os
              << "spline (ratio=" << u << "):" << incindent << iendl
              <<   *spline_ << decindent << iendl
              << "Inner function:" << incindent << iendl
              <<   *inner_ << decindent;
          }

          DifferentiableFunctionPtr_t inner_;
          const value_type u_;

          // LiegroupElement base_;
          // Spline::BasisFunctionVector_t basisFunc;
          const Spline::Ptr_t spline_;
      }; // class SplineFunction

      /// This class computes an output spline from an input spline.
      /// Let \f$ I \f$ be the input spline (0 -> T_i)
      /// and \f$ O \f$ the output spline (0 -> T_o).
      /// - \f$ O(  0) = inner (I(  0)) \f$
      /// - \f$ O(T_o) = inner (I(T_i)) \f$
      /// - \f$ \frac{dO}{dt}(  0) = Jinner (I(  0)) \frac{dI}{dt}(  0) \f$
      /// - \f$ \frac{dO}{dt}(T_o) = Jinner (I(T_i)) \frac{dI}{dt}(T_i) \f$
      template <int _PolynomeBasis, int _SplineOrder>
      class HPP_CORE_LOCAL SplineExplicitFunction :
        public constraints::DifferentiableFunction
      {
        public:
          typedef boost::shared_ptr<SplineExplicitFunction> Ptr_t;
          typedef path::Spline<_PolynomeBasis, _SplineOrder> Spline;

          /// \param base spline base
          /// \param u ratio at which this function applies.
          /// \param T length of the spline. (May be one day this will be an optimization variable.)
          SplineExplicitFunction (const DifferentiableFunctionPtr_t& inner,
              const DevicePtr_t robot,
              const Configuration_t& base,
              const value_type& u,
              const value_type& T) :
            DifferentiableFunction (
                inner->inputDerivativeSize()*Spline::NbCoeffs,
                inner->inputDerivativeSize()*Spline::NbCoeffs,
                inner->outputSpace(), inner->name() + " | " + toStr(inArgs)),
            inner_ (inner),
            // base_ (base, robot->configSpace()->vectorSpaceMerged())
            u_ (u),
            spline_ (Spline::create (robot, interval_t(0, T), ConstraintSetPtr_t()))
          {
            // TODO Compute activeParameters_ and activeDerivativeParameters_
            activeParameters_          .setConstant(true);
            activeDerivativeParameters_.setConstant(true);

            // Spline::basisFunctionDerivative(0, u, basisFunc);
            spline_->base(base_);
          }

        protected:
          typedef Eigen::Map<Spline::ParameterMatrix_t> ParameterMap_t;

          void impl_compute (LiegroupElement& y, vectorIn_t parameters) const
          {
            spline_->rowParameters (parameters);
            Configuration_t q (spline_->outputSize());
            (*spline_) (q, u*spline_->length());
            inner_->value(y, q);
          }

          void impl_jacobian (matrixOut_t J, vectorIn_t parameters) const
          {
            spline_->rowParameters (parameters);
            Configuration_t q (spline_->outputSize());
            const value_type t = u*spline_->length();
            (*spline_) (q, t);

            Eigen::Matrix<value_type, Spline::NbCoeffs, 1> paramDerCoeffs;
            spline_->parameterDerivativeCoefficients (paramDerCoeffs, t);

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
              J.middleCols (i_Pi, Nv).noalias() = paramDerCoeffs(i) * Jinner;
              i_Pi += Nv;
            }
            assert (inputDerivativeSize() == i_Pi);
            // Compute df/dP0
            Jinner *= paramDerCoeffs(0);
          }

          std::ostream& print (std::ostream& os) const
          {
            constraints::DifferentiableFunction::print(os);
            os
              << "spline (ratio=" << u << "):" << incindent << iendl
              <<   *spline_ << decindent << iendl
              << "Inner function:" << incindent << iendl
              <<   *inner_ << decindent;
          }

          DifferentiableFunctionPtr_t inner_;
          const value_type u_;

          // LiegroupElement base_;
          // Spline::BasisFunctionVector_t basisFunc;
          const Spline::Ptr_t spline_;
      }; // class SplineExplicitFunction

      /// The explicit constraint should be handled separately and:
      /// - removed if possible, which means:
      ///   - constant right hand side: the original path treats it explicitely
      ///   - varying right hand side: the initial path treats it explicitely
      /// - made implicit otherwise.

      /// The constraint return is implicit even if an explicit function is
      /// passed.
      /// \param f a function taking as input a robot configuration.
      /// TODO this function should take 
      template <int _PB, int _SO>
      constraints::ImplicitPtr_t applyConstraintOnSpline (
          const DevicePtr_t& device,
          const constraints::ImplicitPtr_t& ic,
          const std::vector<Spline<_PB, _SO>::Ptr_t>& splines,
          const size_type& iSpline,
          const value_type & u)
      {
        assert (f->inputSize          () == device->configSize());
        assert (f->inputDerivativeSize() == device->numerDof());

        typedef SplineFunction<_PB, _SO> SplineFunction_t;
        typedef SplineFunction_t::Spline Spline;
        using function::OfParameterSubset;
        using function::OfParameterSubsetPtr_t;

        SplinePtr_t spline (splines[iSpline]);
        SplineFunction::Ptr sf (new SplineFunction (ic->function(),
              device, spline->base(), u, spline->length()));

        const size_type nDofSpline = Spline::NbCoeffs * device->numberDof(),
              OfParameterSubsetPtr_t ops (OfParameterSubset::create (sf,
                    splines.size() * nDofSpline,
                    splines.size() * nDofSpline,
                    segment_t(iSpline * nDofSpline, nDofSpline),
                    segment_t(iSpline * nDofSpline, nDofSpline)));

        constraints::ExplicitPtr_t ec = HPP_DYNAMIC_PTR_CAST(
            constraints::Explicit, ic);
        return constraints::Implicit::create (ops, ic->comparisonTypes(),
            ic->rightHandSide());
      }
    } // namespace pathOptimization
  } // namespace core
} // namespace hpp
