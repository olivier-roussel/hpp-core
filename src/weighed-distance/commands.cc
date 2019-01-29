// Copyright (c) 2019, Joseph Mirabel
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

#include <../src/weighed-distance/commands.hh>

#include <hpp/core/problem.hh>
#include <hpp/core/problem-solver.hh>

namespace hpp {
  namespace core {
    namespace distance {
      Parameter GetWeights::doExecute(ProblemSolverPtr_t ps)
      {
        if (!ps->problem()) throw std::runtime_error ("The problem is not defined.");
        WeighedDistancePtr_t dist = HPP_DYNAMIC_PTR_CAST(WeighedDistance,
            ps->problem()->distance());
        if (!dist) throw std::runtime_error ("Distance is not of type WeighedDistance.");

        return Parameter (dist->weights ());
      }

      Parameter SetWeights::doExecute(ProblemSolverPtr_t ps)
      {
        if (!ps->problem()) throw std::runtime_error ("The problem is not defined.");
        WeighedDistancePtr_t dist = HPP_DYNAMIC_PTR_CAST(WeighedDistance,
            ps->problem()->distance());
        if (!dist) throw std::runtime_error ("Distance is not of type WeighedDistance.");

        vector_t ws = getParameterValues()[0].vectorValue ();
        dist->weights (ws);

        return Parameter ();
      }
    } // namespace distance
  } // namespace core
} // namespace hpp
