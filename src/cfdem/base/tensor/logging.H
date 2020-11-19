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

#include <inttypes.h>
#include <iostream>
#include <cstdlib> // abort()

namespace base {

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

/*! \brief double type that will be used in default by cfdem */
typedef double d_real_t;

/*! \brief device name CPU */
struct cpu {
  /*! \brief whether this device is CPU or not */
  static const bool kDevCPU = true;
  /*! \brief device flag number, identifies this device */
  static const int kDevMask = 1 << 0;
};

/*! \brief device name GPU */
struct gpu {
  /*! \brief whether this device is CPU or not */
  static const bool kDevCPU = false;
  /*! \brief device flag number, identifies this device */
  static const int kDevMask = 1 << 1;
};

} // namespace base

/*! \brief type that will be used for index */
#ifndef CFDEM_INT64_TENSOR_SIZE
#define CFDEM_INT64_TENSOR_SIZE 0
#endif

/*! \brief whether do padding during allocation */
#ifndef CFDEM_ALLOC_PAD
#define CFDEM_ALLOC_PAD true
#endif

/*! \brief whether use SSE */
#ifndef CFDEM_USE_SSE
  #define CFDEM_USE_SSE 1
#endif

#if CFDEM_INT64_TENSOR_SIZE == 1
  typedef int64_t index_t;
#else
  typedef int32_t index_t;
#endif

/*! \brief default data type for tensor string */
#ifndef CFDEM_DEFAULT_DTYPE
#define CFDEM_DEFAULT_DTYPE = ::base::d_real_t
#endif

/*! \brief default device type */
#ifndef CFDEM_DEFAULT_DEVICE
#define CFDEM_DEFAULT_DEVICE = ::base::cpu
#endif

#define CHECK(x)                                                                                                   \
  if (!(x))                                                                                                        \
    base::LogMessageFatal(__FILE__, __LINE__).stream() << "Check (" << #x << ") faild! "
#define CHECK_GT(x, y) CHECK((x) > (y))
#define CHECK_LE(x, y) CHECK((x) <= (y))
#define CHECK_GE(x, y) CHECK((x) >= (y))
#define CHECK_EQ(x, y) CHECK((x) == (y))
#define CHECK_NE(x, y) CHECK((x) != (y))
#define CHECK_NOTNULL(x)                                                                                           \
  ((x) == NULL ? base::LogMessageFatal(__FILE__, __LINE__).stream() << "Check not null failed: " #x, (x) : (x))

#endif // __LOGGING_H__
