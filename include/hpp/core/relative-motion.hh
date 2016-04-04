//
// Copyright (c) 2016 CNRS
// Authors: Joseph Mirabel
//
// This file is part of hpp-core
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
// hpp-core  If not, see
// <http://www.gnu.org/licenses/>.

#ifndef HPP_CORE_RELATIVE_MOTION_HH
#define HPP_CORE_RELATIVE_MOTION_HH

#include <Eigen/Core>

#include <hpp/model/fwd.hh>
#include <hpp/core/fwd.hh>

namespace hpp {
  namespace core {
    struct RelativeMotion {
      enum RelativeMotionType {
        /// The relative motion is fully constrained and the constraint cannot be
        /// parameterized (has constant right-hand side)
        Constrained,
        /// The relative motion is fully constrained but the constraint can be
        /// parameterized (has non-constant right-hand side)
        Parameterized,
        /// The relative motion is not constrained
        Unconstrained
      };

      typedef Eigen::Matrix<RelativeMotionType, Eigen::Dynamic, Eigen::Dynamic> matrix_type;

      static matrix_type matrix (const DevicePtr_t& robot);

      /// \todo LockedJoint always has a non-constant RHS which means it will
      ///       always be treated a parameterized constraint. Even when the
      ///       value is not going to change...
      static void fromConstraint (
          matrix_type& matrix,
          const DevicePtr_t& robot,
          const ConstraintSetPtr_t& constraint);
    };
  } // namespace core
} // namespace hpp

namespace Eigen {
  template<> struct NumTraits<hpp::core::RelativeMotion::RelativeMotionType>
    : NumTraits<int> {};
} // namespace Eigen

#endif // HPP_CORE_RELATIVE_MOTION_HH
