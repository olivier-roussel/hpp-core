//
// Copyright 2010-2018 CNRS
//
// Author: Florent Lamiraux, Joseph Mirabel
//
// This file is part of hpp-core.
// hpp-core is free software: you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
// hpp-core is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.  You should
// have received a copy of the GNU Lesser General Public License along
// with hpp-core.  If not, see <http://www.gnu.org/licenses/>.

#include <sstream>

#include <hpp/util/exception-factory.hh>

#include <hpp/core/command.hh>
#include <hpp/core/problem-solver.hh>

namespace hpp {
  namespace core {

    Command::Command(const std::string& docstring) :
      valueTypeVector_(), docstring_(docstring)
    {
    }

    void Command::valueTypes (const std::vector<Parameter::Type>& valueTypes)
    {
      valueTypeVector_ = valueTypes;
    }

    const std::vector<Parameter::Type>& Command::valueTypes() const
    {
      return valueTypeVector_;
    }

    void Command::setParameterValues(const std::vector<Parameter>& values)
    {
      const std::vector<Parameter::Type>& paramTypes = valueTypes();
      // Check that number of parameters is correct
      if (values.size() != paramTypes.size()) {
        throw std::invalid_argument ("wrong number of parameters");
      }
      // Check that each parameter is of correct type
      for (unsigned int iParam=0; iParam < values.size(); iParam++) {
        if (values[iParam].type() != paramTypes[iParam]) {
          HPP_THROW (std::invalid_argument, "Argument " << iParam <<
              " is of wrong type: " << Parameter::typeName(paramTypes[iParam])
              << " expected, got " << Parameter::typeName(values[iParam].type())
              );
        }
      }
      // Copy vector of values in private part
      valueVector_ = values;
    }

    const std::vector<Parameter>& Command::getParameterValues() const
    {
      return valueVector_;
    }

    Parameter Command::execute(ProblemSolverPtr_t ps)
    {
      return doExecute(ps);
    }

    std::string Command::getDocstring() const
    {
      return docstring_;
    }
  } // namespace core
} //namespace hpp
