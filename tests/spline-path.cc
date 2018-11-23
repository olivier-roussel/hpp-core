// Copyright (c) 2017, Joseph Mirabel
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

#define BOOST_TEST_MODULE spline_path
#include <boost/test/included/unit_test.hpp>

#include <hpp/core/path/spline.hh>
#include <hpp/core/path/composite-bezier.hh>

#include <pinocchio/algorithm/joint-configuration.hpp>

#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/joint-collection.hh>
#include <hpp/pinocchio/simple-device.hh>
#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/liegroup.hh>

#include <hpp/core/problem.hh>
#include <hpp/core/steering-method/straight.hh>
#include <hpp/core/steering-method/spline.hh>

#include <../tests/util.hh>

#define TOSTR( x ) static_cast< std::ostringstream & >( ( std::ostringstream() << x ) ).str()

using namespace hpp::core;
using namespace hpp::pinocchio;

DevicePtr_t createRobot ()
{
  DevicePtr_t robot = unittest::makeDevice(unittest::HumanoidRomeo);
  robot->controlComputation((Computation_t) (JOINT_POSITION | JACOBIAN));
  robot->rootJoint()->lowerBound (0, -1);
  robot->rootJoint()->lowerBound (1, -1);
  robot->rootJoint()->lowerBound (2, -1);
  robot->rootJoint()->upperBound (0,  1);
  robot->rootJoint()->upperBound (1,  1);
  robot->rootJoint()->upperBound (2,  1);
  return robot;
}

typedef std::pair<value_type, value_type> Pair_t;

std::ostream& operator<< (std::ostream& os, const Pair_t& p) {
  os << "Pair " << p.first << ", " << p.second;
  return os;
}

void printAt(const PathPtr_t& p, ConfigurationOut_t q, value_type t)
{
  (*p)(q, t);
  std::cout << t << ":\t" << q.transpose() << std::endl;
}
void checkAt(const PathPtr_t orig, value_type to,
             const PathPtr_t extr, value_type te) {
  Configuration_t q1 (orig->outputSize()),
                  q2 (orig->outputSize());
  (*orig)(q1, to);
  (*extr)(q2, te);
  BOOST_CHECK_MESSAGE(q2.isApprox(q1),
      "\nPath 1: " << q1.head<10>().transpose() <<
      "\nPath 2: " << q2.head<10>().transpose()
      );
}

template <int SplineType> void compare_to_straight_path ()
{
  typedef path::Spline<SplineType, 1> path_t;
  typedef steeringMethod::Spline<SplineType, 1> SM_t;

  DevicePtr_t dev = createRobot();
  BOOST_REQUIRE (dev);
  Problem problem (dev);

  Configuration_t q1 (::pinocchio::randomConfiguration(dev->model()));
  Configuration_t q2 (::pinocchio::randomConfiguration(dev->model()));

  vector_t v (dev->numberDof());
  difference<RnxSOnLieGroupMap> (dev, q2, q1, v);

  // create StraightPath
  PathPtr_t sp = (*problem.steeringMethod()) (q1, q2);
  // value_type length = sp->length();

  // Create linear spline
  typename SM_t::Ptr_t sm (SM_t::create (problem));
  PathPtr_t ls_abstract = (*sm) (q1, q2);
  typename path_t::Ptr_t ls = HPP_DYNAMIC_PTR_CAST(path_t, ls_abstract);

  /*
  typename path_t::Ptr_t ls = path_t::create (dev, interval_t(0, length), ConstraintSetPtr_t());
  ls->base (q1);
  typename path_t::ParameterMatrix_t ls_param = ls->parameters();
  ls_param.row(0).setZero();
  ls_param.row(1) = v;
  ls->parameters(ls_param);
  */

  CONFIGURATION_VECTOR_IS_APPROX(dev, sp->initial(), ls->initial(), 1e-7);
  CONFIGURATION_VECTOR_IS_APPROX(dev, sp->end()    , ls->end()    , 1e-7);

  const size_type N = 10;
  const value_type step1 = sp->length() / N;
  const value_type step2 = ls->length() / N;

  // Check that straight path and linear spline return the same result.
  for (size_type i = 0; i < N; ++i)
    checkAt (sp, value_type(i) * step1, ls, value_type(i) * step2);

  // Check that the velocities are equals
  vector_t v1 (dev->numberDof());
  vector_t v2 (dev->numberDof());
  for (size_type i = 0; i < N; ++i) {
    sp->derivative(v1, value_type(i) * step1, 1);
    ls->derivative(v2, value_type(i) * step2, 1);
    BOOST_CHECK_SMALL((v2 * step2 - v1 * step1).squaredNorm(), 1e-12);
  }

  // Check integral value
  BOOST_CHECK_EQUAL(ls->squaredNormIntegral(2), 0);
  BOOST_CHECK_CLOSE(ls->squaredNormIntegral(1), v.squaredNorm() / sp->length(), 1e-12);

  vector_t derivative (dev->numberDof() * 2);
  ls->squaredNormIntegralDerivative(1, derivative);
  switch (SplineType) {
    case path::CanonicalPolynomeBasis:
      BOOST_CHECK( derivative.head(dev->numberDof()).isZero());
      EIGEN_VECTOR_IS_APPROX(derivative.tail(dev->numberDof()), (2./sp->length()) * v , 1e-6);
      break;
    case path::BernsteinBasis:
      EIGEN_VECTOR_IS_APPROX(derivative.head(dev->numberDof()), (-2./sp->length()) * v, 1e-6);
      EIGEN_VECTOR_IS_APPROX(derivative.tail(dev->numberDof()), ( 2./sp->length()) * v, 1e-6);
      break;
  } 
}

template <int SplineType, int Degree>
void check_velocity_bounds ()
{
  typedef steeringMethod::Spline<SplineType, Degree> SM_t;

  DevicePtr_t dev = createRobot();
  BOOST_REQUIRE (dev);
  Problem problem (dev);

  Configuration_t q1 (::pinocchio::randomConfiguration(dev->model()));
  Configuration_t q2 (::pinocchio::randomConfiguration(dev->model()));
  std::vector<int> orders (1, 1);
  vector_t v1 (vector_t::Random(dev->numberDof())),
           v2 (vector_t::Random(dev->numberDof()));


  // Create spline
  typename SM_t::Ptr_t sm (SM_t::create (problem));
  PathPtr_t spline = sm->steer (q1, orders, v1, q2, orders, v2);

  vector_t vb1 (vector_t::Random(dev->numberDof())), vb2 = vb1;
  value_type t0 = spline->timeRange().first, t1 = spline->timeRange().second;
  spline->velocityBound (vb1, t0, t1);

  std::size_t N = 1000;
  value_type step = spline->length() / value_type(N);
  for (std::size_t i = 0; i < N; ++i) {
    spline->velocityBound (vb2, t0, t1);
    BOOST_CHECK_MESSAGE((vb2.array() <= vb1.array()).all(),
        "i=" << i << " Velocity bound should have decreased. Interval is ["
        << t0 << ", " << t1 << "]. Difference:\n"
        << (vb1 - vb2).transpose());
    vb1 = vb2;
    if (i%2) {
      t0 += step;
    } else {
      t1 -= step;
    }
  }
}

BOOST_AUTO_TEST_CASE (spline_canonical)
{
  compare_to_straight_path<path::CanonicalPolynomeBasis>();
}

BOOST_AUTO_TEST_CASE (spline_bernstein)
{
  compare_to_straight_path<path::BernsteinBasis>();
}

BOOST_AUTO_TEST_CASE (spline_bernstein_velocity)
{
  check_velocity_bounds<path::BernsteinBasis, 3>();
}

BOOST_AUTO_TEST_CASE (composite_bezier_1)
{
  typedef path::CompositeBezier<1> CompositeBezier;
  typedef typename CompositeBezier::Ptr_t CompositeBezierPtr_t;

  const size_type N = 1;
  LiegroupSpacePtr_t space = LiegroupSpace::Rn(N);
  Configuration_t q1 (N), q2 (N);
  vector_t v (N);
  q1 << 0;
  q2 << 1;
  CompositeBezierPtr_t path = CompositeBezier::create (space, q1, q2, 1);

  (*path) (q1, 0.5);
  BOOST_CHECK_EQUAL(q1[0], 0.5);

  path->derivative (v, 0.5, 1);
  BOOST_CHECK_EQUAL(v[0], 1);
}

BOOST_AUTO_TEST_CASE (composite_bezier_3)
{
  typedef path::CompositeBezier<3> CompositeBezier;
  typedef typename CompositeBezier::Ptr_t CompositeBezierPtr_t;

  const size_type N = 1;
  LiegroupSpacePtr_t space = LiegroupSpace::Rn(N);
  Configuration_t q1 (N), q2 (N);
  vector_t v1 (N), v2 (N);
  q1 << 0;
  q2 << 1;
  // Path (t) = t*t*(3-2*t)
  // Path' (t) = 6*t*(1-t)
  CompositeBezierPtr_t path = CompositeBezier::create (space, q1, q2, 1);

  matrix_t controlPoints (N, 4);
  controlPoints << 0, 0, 1, 1;
  BOOST_CHECK_EQUAL (controlPoints, path->controlPoints());

  (*path) (q1, 0.5);
  BOOST_CHECK_EQUAL(q1[0], 0.5);

  value_type t = 0;
  path->derivative (v1, t, 1);
  BOOST_CHECK_EQUAL(v1[0], 6*t*(1-t));

  t = 0.5;
  path->derivative (v1, t, 1);
  BOOST_CHECK_EQUAL(v1[0], 6*t*(1-t));

  t = 1;
  path->derivative (v1, t, 1);
  BOOST_CHECK_EQUAL(v1[0], 6*t*(1-t));

  v1 << 0.25;
  path->velocityAtControlPoint(0, v1);
  path->velocityAtControlPoint(1, v1);

  const value_type tol = 1e-4;
  path->derivative (v2, 0, 1);
  BOOST_CHECK_CLOSE(v1[0], v2[0], tol);
  path->derivative (v2, 1, 1);
  BOOST_CHECK_CLOSE(v1[0], v2[0], tol);

  q1 << -0.5;
  std::cout << path->controlPoints() << std::endl;
  const std::vector<value_type>& params = path->params();
  for (std::size_t i = 0; i < params.size(); ++i)
    std::cout << params[i] << ", ";
  std::cout << std::endl;

  size_type idx = path->insert (0.5, q1);
  BOOST_CHECK_EQUAL(idx, 1);
  std::cout << path->controlPoints() << std::endl;
  //std::copy(path->params().begin(), path->params().end(), std::ostream_iterator<char>(std::cout, ", "));
  for (std::size_t i = 0; i < params.size(); ++i)
    std::cout << params[i] << ", ";
  std::cout << std::endl;

  path->velocityAtControlPoint(0, v1);
  path->velocityAtControlPoint(2, v1);

  path->derivative (v2, 0, 1);
  BOOST_CHECK_CLOSE(v1[0], v2[0], tol);
  path->derivative (v2, 1, 1);
  BOOST_CHECK_CLOSE(v1[0], v2[0], tol);
}
