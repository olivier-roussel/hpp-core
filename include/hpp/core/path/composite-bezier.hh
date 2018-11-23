// Copyright (c) 2018 CNRS
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

#ifndef HPP_CORE_PATH_COMPOSITE_BEZIER_HH
# define HPP_CORE_PATH_COMPOSITE_BEZIER_HH

# include <boost/core/enable_if.hpp>

# include <hpp/pinocchio/liegroup-element.hh>

# include <hpp/core/fwd.hh>
# include <hpp/core/config.hh>
# include <hpp/core/path.hh>

namespace hpp {
  namespace core {
    namespace path {
      /// \addtogroup path
      /// \{

      /// Piecewise Bezier curve.
      ///
      /// This type of path is the return type of PathProjector algorithms.
      ///
      /// Degrees of freedom are interpolated depending on the type of
      /// \link hpp::pinocchio::Joint joint \endlink
      /// they parameterize:
      ///   \li linear interpolation for translation joints, bounded rotation
      ///       joints, and translation part of freeflyer joints,
      ///   \li angular interpolation for unbounded rotation joints,
      ///   \li constant angular velocity for SO(3) part of freeflyer joints.
      template <int Order>
      class HPP_CORE_DLLAPI CompositeBezier : public Path
      {
      public:
        typedef boost::shared_ptr<CompositeBezier> Ptr_t;

        typedef std::vector<value_type> Params_t;
        /// Each column corresponds to one control point.
        typedef matrix_t ControlPoints_t;
        typedef matrix_t::ConstColXpr ControlPointConst_t;
        typedef matrix_t::ColXpr      ControlPoint_t;

        /// Destructor
        virtual ~CompositeBezier () throw () {}

        /// Create instance and return shared pointer
        /// \param space output Lie Group space
        /// \param init, end Start and end configurations of the path
        /// \param length Distance between the configurations.
        static Ptr_t create (const LiegroupSpacePtr_t& space,
                                         ConfigurationIn_t init,
                                         ConfigurationIn_t end,
                                         value_type length)
        {
          CompositeBezier* ptr = new CompositeBezier (space, init, end, length);
          Ptr_t shPtr (ptr);
          ptr->init (shPtr);
          return shPtr;
        }

        /// Create instance and return shared pointer
        /// \param init, end Start and end configurations of the path
        /// \param length Distance between the configurations.
        static Ptr_t create (pinocchio::LiegroupElementConstRef init,
                             pinocchio::LiegroupElementConstRef end,
                             value_type length)
        {
          assert (init.space() == end.space());
          CompositeBezier* ptr = new CompositeBezier (init.space(),
              init.vector(), end.vector(), length);
          Ptr_t shPtr (ptr);
          ptr->init (shPtr);
          return shPtr;
        }

        /// Create instance and return shared pointer
        /// \param space Robot corresponding to configurations
        /// \param init, end Start and end configurations of the path
        /// \param length Distance between the configurations.
        /// \param constraints the path is subject to
        static Ptr_t create (const LiegroupSpacePtr_t& space,
                                         ConfigurationIn_t init,
                                         ConfigurationIn_t end,
                                         value_type length,
                                         ConstraintSetPtr_t constraints)
        {
          CompositeBezier* ptr = new CompositeBezier (space, init, end, length,
                                                constraints);
          Ptr_t shPtr (ptr);
          ptr->init (shPtr);
          return shPtr;
        }

        /// Return a shared pointer to a copy of this
        virtual PathPtr_t copy () const
        {
          CompositeBezier* ptr = new CompositeBezier (*this);
          Ptr_t shPtr (ptr);
          ptr->initCopy (shPtr);
          return shPtr;
        }

        /// Return a shared pointer to a copy of this and set constraints
        ///
        /// \param constraints constraints to apply to the copy
        /// \precond *this should not have constraints.
        virtual PathPtr_t copy (const ConstraintSetPtr_t& constraints) const
        {
          CompositeBezier* ptr = new CompositeBezier (*this, constraints);
          Ptr_t shPtr (ptr);
          ptr->initCopy (shPtr);
          return shPtr;
        }

        virtual PathPtr_t reverse () const;

        /// Return the internal robot.
        LiegroupSpacePtr_t space () const
        {
          return space_;
        }

        /// Insert a control point
        /// by adding \f$ 2 * Order + 1 \f$ points at the corresponding place.
        /// \return the index of the control point.
        size_type insert (const value_type& param, ConfigurationIn_t config);

        /// Set velocity
        /// Only valid for Order >= 3.
        void velocityAtControlPoint (const size_type& idx, vectorIn_t v);

        /// Get the initial configuration
        Configuration_t initial () const
        {
          return controlPoints_.col(0);
        }

        /// Get the final configuration
        Configuration_t end () const
        {
          return controlPoints_.col(npts_-1);
        }

        ControlPoints_t controlPoints () const
        {
          return controlPoints_.leftCols(npts_);
        }

        const Params_t& params () const
        {
          return params_;
        }

      protected:
        /// Print path in a stream
        virtual std::ostream& print (std::ostream &os) const
        {
          os << "CompositeBezier:" << std::endl;
          Path::print (os);
          os << "initial configuration: " << initial().transpose () << std::endl;
          os << "final configuration:   " << end().transpose () << std::endl;
          return os;
        }

        /// Constructor
        CompositeBezier (const LiegroupSpacePtr_t& space, ConfigurationIn_t init,
                      ConfigurationIn_t end, value_type length);

        /// Constructor with constraints
        CompositeBezier (const LiegroupSpacePtr_t& space, ConfigurationIn_t init,
                      ConfigurationIn_t end, value_type length,
                      ConstraintSetPtr_t constraints);

        /// Copy constructor
        CompositeBezier (const CompositeBezier& path);

        /// Copy constructor with constraints
        CompositeBezier (const CompositeBezier& path,
                      const ConstraintSetPtr_t& constraints);

        void init (Ptr_t self);

        void initCopy (Ptr_t self);

        virtual bool impl_compute (ConfigurationOut_t result,
                                   value_type param) const;
        /// Virtual implementation of derivative
        virtual void impl_derivative (vectorOut_t result, const value_type& t,
                                      size_type order) const;
        virtual void impl_velocityBound (vectorOut_t result,
            const value_type& t0, const value_type& t1) const;

        /// Extraction/Reversion of a sub-path
        /// See Path::extract
        //PathPtr_t impl_extract (const interval_t& subInterval) const
          //throw (projection_error);
      private:
        typedef boost::weak_ptr<CompositeBezier> WkPtr_t;

        /// Constructor from control points
        template <typename Derived>
        CompositeBezier (const LiegroupSpacePtr_t& space,
            const Params_t& params,
            const Eigen::DenseBase<Derived>& controlPoints,
            ConstraintSetPtr_t constraints);

        size_type indexAtParam (const value_type& param,
          value_type& T, value_type& localParam) const;

        LiegroupSpacePtr_t space_;
        Params_t params_;
        ControlPoints_t controlPoints_;
        size_type npts_;
        WkPtr_t weak_;
      }; // class CompositeBezier
      /// \}
    } //   namespace path
  } //   namespace core
} // namespace hpp
#endif // HPP_CORE_PATH_COMPOSITE_BEZIER_HH
