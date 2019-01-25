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

#ifndef HPP_CORE_COMMAND_HH
#define HPP_CORE_COMMAND_HH

#include <vector>
#include <hpp/core/fwd.hh>
#include <hpp/core/config.hh>
#include <hpp/core/parameter.hh>

#include <boost/assign/list_of.hpp>

namespace hpp {
  namespace core {
    using boost::assign::list_of;

    /// Abstract class for commands stored in ProblemSolver
    ///
    /// This class provide a mean to control entities from external python script.
    ///
    /// A command
    /// \li is owned by the ProblemSolver,
    /// \li takes parameters of type Parameter,
    /// \li return an instance of Parameter when calling Command::execute()
    ///
    /// At construction, the prototype of the command is defined by providing
    /// a vector of Parameter::Type.
    ///
    /// Parameters are set by calling Command::setParameterValues with a
    /// vector of Values the types of which should fit the vector specified
    /// at construction.
    class HPP_CORE_DLLAPI Command
    {
    public:
      virtual ~Command() {}
      /// Creates an command without arguments
      /// \param docstring documentation of the command
      Command(const std::string& docstring);
      /// Store a vector of value types
      /// \param valueTypes vector specifying the number and types of parameters
      void valueTypes (const std::vector<Parameter::Type>& valueTypes);
      /// Return the value type of all parameters
      const std::vector<Parameter::Type>& valueTypes() const;
      /// Set parameter values
      void setParameterValues(const std::vector<Parameter>& values);
      /// Get parameter values
      const std::vector<Parameter>& getParameterValues() const;
      /// Execute the command after checking parameters
      Parameter execute(ProblemSolverPtr_t ps);
      /// Get documentation string
      std::string getDocstring() const;
    protected:
      /// Specific action performed by the command
      virtual Parameter doExecute(ProblemSolverPtr_t ps) = 0;
    private:
      std::vector<Parameter::Type> valueTypeVector_;
      std::vector<Parameter> valueVector_;
      std::string docstring_;
    };
  } // namespace core
} // namespace hpp

#endif //HPP_CORE_COMMAND_HH
