/*!
 * Copyright (c) 2020 by Contributors
 * @file int32.C
 *
 * @author Wang Ding
 */

#include "./int32.H"
#include "IOstreams.H"

#include <inttypes.h>
#include <sstream>
#include <cerrno>

namespace base {

std::string toString(const int32_t val) {
  std::ostringstream oss;
  oss << val;
  return oss.str();
}

int32_t readInt32(Istream& is) {
    int32_t val;
    is >> val;
    return val;
}

} // namespace base
