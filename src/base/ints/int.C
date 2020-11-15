/*!
 * Copyright (c) 2020 by Contributors
 * @file int.C
 * \brief system int
 *
 * @author Wang Ding
 */

#ifndef __INT_H__
#define __INT_H__

#include "int.H"
#include "IOstreams.H"

namespace base {

int readInt(Istream& is) {
  int val;
  is >> val;
  return val;
}

} // namespace base

#endif // __INT_H__
