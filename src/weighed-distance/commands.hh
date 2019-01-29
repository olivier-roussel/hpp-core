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

#ifndef HPP_CORE_WEIGHED_DISTANCE_COMMANDS_HH
#define HPP_CORE_WEIGHED_DISTANCE_COMMANDS_HH

#include <hpp/core/config.hh>
#include <hpp/core/command.hh>
#include <hpp/core/weighed-distance.hh>

namespace hpp {
  namespace core {
    namespace distance {
      class HPP_CORE_LOCAL GetWeights : public Command
      {
        public:
          GetWeights ()
            : Command ("Get weights of a WeighedDistance")
          {}

        protected:
          Parameter doExecute(ProblemSolverPtr_t ps);
      };
      class HPP_CORE_LOCAL SetWeights : public Command
      {
        public:
          SetWeights ()
            : Command ("Set weights of a WeighedDistance")
          {
            valueTypes (list_of(Parameter::VECTOR));
          }

        protected:
          Parameter doExecute(ProblemSolverPtr_t ps);
      };
    } // namespace distance
  } // namespace core
} // namespace hpp

#endif // HPP_CORE_WEIGHED_DISTANCE_COMMANDS_HH
