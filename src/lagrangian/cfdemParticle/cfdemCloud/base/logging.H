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
    Defines logging of cfdem

File
    logging.H
\*---------------------------------------------------------------------------*/

#ifndef __LOGGING_H__
#define __LOGGING_H__

#include <iostream>

namespace cfdem {

class LogMessage {
public:
  LogMessage(const char* file, int line): log_stream_(std::cerr) {
    log_stream_ << file << ": " << line << ": ";
  }
  ~LogMessage() {
    log_stream_ << std::endl;
  }
  inline std::ostream& stream() { return log_stream_; }

protected:
  std::ostream& log_stream_;

private:
  LogMessage(const LogMessage&);
  LogMessage& operator= (const LogMessage&);
};

class LogMessageFatal: public LogMessage {
public:
  LogMessageFatal(const char* file, int line): LogMessage(file, line) {}
  ~LogMessageFatal() {
    abort();
  }

private:
  LogMessageFatal(const LogMessageFatal&);
  LogMessageFatal& operator= (const LogMessageFatal&);
};

#define CHECK(x) \
  if (!(x)) { \
    cfdem::LogMessageFatal(__FILE__, __LINE__).stream() << "Check failed: " #x; \
  }
#define CHECK_LT(x, y) CHECK((x) < (y))
#define CHECK_GT(x, y) CHECK((x) > (y))
#define CHECK_LE(x, y) CHECK((x) <= (y))
#define CHECK_GE(x, y) CHECK((x) >= (y))
#define CHECK_EQ(x, y) CHECK((x) == (y))
#define CHECK_NE(x, y) CHECK((x) != (y))
#define CHECK_NOTNULL(x) \
  ((x) == NULL ? cfdem::LogMessageFatal(__FILE__, __LINE__).stream() << "Check not null failed: " #x, (x) : (x))
} // namespace cfdem

#endif // __LOGGING_H__
