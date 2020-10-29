/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Basic functions and classes for exceptions.

File
    exceptions.H
\*---------------------------------------------------------------------------*/

#ifndef __EXCEPTIONS_H__
#define __EXCEPTIONS_H__

#include <stdexcept> // for std::runtime_error
#include <vector>
#include "logging.H"

namespace cfdem {

enum ExceptionType {
  kTypeNull = -1,
  kIndexOutOfRange = 0,
  kIllegalParticleRatio = 1
};

static std::vector<std::string> ExceptionMessages = {
  "Index out of range: ",
  "Illegal particle ratio: "
};

class RuntimeException : public std::runtime_error {
public:
  explicit RuntimeException(const std::string& message,
                            ExceptionType ex_type = ExceptionType::kTypeNull)
    : ex_type_(ex_type), std::runtime_error(message) {}

  const char* what() const noexcept(true) override {
    if (ex_type_ < 0) {
      return std::runtime_error::what();
    }
    CHECK_LT(ex_type_, ExceptionMessages.size())
    std::string message(std::runtime_error::what());
    return (ExceptionMessages[ex_type_] + message).c_str();
  }

private:
  ExceptionType ex_type_;
};

} // namespace cfdem

/*!
 * \brief Run function and catch error, log unknown error.
 * \param func Expression to call.
 */
#define CFDEM_CATCH_ERROR(func)                                      \
{                                                                    \
  try {                                                              \
    (func);                                                          \
  } catch(const RuntimeException& e) {                               \
    Pout << e.what() << endl;                                        \
  }                                                                  \
}

/*!
 * \brief Run function and catch error, log unknown error.
 * \param func Expression to call.
 */
#define CFDEM_CATCH_ERROR_RETURN(func)                               \
{                                                                    \
  try {                                                              \
    return (func);                                                   \
  } catch(const RuntimeException& e) {                               \
    Pout << e.what() << endl;                                        \
  }                                                                  \
}

#endif // __EXCEPTIONS_H__
