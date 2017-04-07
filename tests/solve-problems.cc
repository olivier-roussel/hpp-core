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

#define BOOST_TEST_MODULE solve_problems
#include <boost/test/included/unit_test.hpp>

// Force benchmark output
#define HPP_ENABLE_BENCHMARK 1
#include <hpp/util/timer.hh>

#include <hpp/pinocchio/simple-device.hh>
#include <hpp/pinocchio/device.hh>

#include <hpp/core/path-vector.hh>
#include <hpp/core/problem-solver.hh>

using namespace hpp::core;

HPP_DEFINE_TIMECOUNTER(solve);
HPP_DEFINE_TIMECOUNTER(optimize);

BOOST_AUTO_TEST_CASE (ur_benchmark)
{
  ProblemSolverPtr_t ps = ProblemSolver::create();

  DevicePtr_t robot = hpp::pinocchio::unittest::makeDevice(hpp::pinocchio::unittest::ManipulatorArm2);
  ps->robot(robot);
  ConfigurationPtr_t qinit (new Configuration_t(robot->neutralConfiguration()));
  ConfigurationPtr_t qgoal (new Configuration_t(*qinit));
  *qgoal << 0.019000000000000003, -1.4470000000000003, 1.047, -2.809845514297485, 0.9638400111347437, 0.0, 1.2144489479064942, 0.0, 0.02, -0.02, 1.013, 0.0, 0.427585186958313, 0.7770800089463592, -0.42826000213623044, 1.8374643480777741, 0.0, 0.02, -0.02; 

  /*
  ps->pathValidationType("Discretized", 0.05);
  ps->pathPlannerType("DiffusingPlanner");

  ps->initConfig(qinit);
  ps->addGoalConfig(qgoal);

  HPP_START_TIMECOUNTER (solve);
  ps->solve();
  HPP_STOP_TIMECOUNTER (solve);
  HPP_STREAM_TIMECOUNTER (std::cout, solve) << std::endl;
  std::cout << "Solution length " << ps->paths()[0]->length() << std::endl;

  ps->addPathOptimizer("RandomShortcut");
  HPP_START_TIMECOUNTER (optimize);
  ps->optimizePath(ps->paths()[0]);
  HPP_STOP_TIMECOUNTER (optimize);
  HPP_STREAM_TIMECOUNTER (std::cout, optimize) << std::endl;

  std::cout << "Optimized length " << ps->paths()[1]->length() << std::endl;
*/
  ps = ProblemSolver::create();
  ps->robot(robot);
  ps->pathValidationType("Dichotomy", 0);
  ps->pathPlannerType("DiffusingPlanner");

  ps->initConfig(qinit);
  ps->addGoalConfig(qgoal);

  HPP_START_TIMECOUNTER (solve);
  ps->solve();
  HPP_STOP_TIMECOUNTER (solve);
  HPP_STREAM_TIMECOUNTER (std::cout, solve) << std::endl;
  std::cout << "Solution length " << ps->paths()[0]->length() << std::endl;

  ps->addPathOptimizer("GradientBased");
  HPP_START_TIMECOUNTER (optimize);
  ps->optimizePath(ps->paths()[0]);
  HPP_STOP_TIMECOUNTER (optimize);
  HPP_STREAM_TIMECOUNTER (std::cout, optimize) << std::endl;

  std::cout << "Optimized length " << ps->paths()[1]->length() << std::endl;
}
