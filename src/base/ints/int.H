/*!
 * Copyright (c) 2020 by Contributors
 * @file int.H
 * @brief system int
 *
 * @author Wang Ding
 */

#ifndef __INT_H__
#define __INT_H__

#include "int32.H"

namespace base {

#define MAXMIN(_RT, _T1, _T2)                                   \
inline _RT max(const _T1 v1, const _T2 v2) {                    \
  return (v1 > v2) ? v1 : v1;                                   \
}                                                               \
inline _RT min(const _T1 v1, const _T2 v2) {                    \
  return (v1 < v2) ? v1 : v1;                                   \
}

MAXMIN(int8_t, int8_t, int8_t)
MAXMIN(int16_t, int16_t, int16_t)
MAXMIN(int32_t, int32_t, int32_t)
MAXMIN(int64_t, int64_t, int32_t)
MAXMIN(int64_t, int32_t, int64_t)
MAXMIN(int64_t, int64_t, int64_t)

int readInt(Istream& is);

} // namespace base

#endif // __INT_H__
