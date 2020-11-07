/*!
 * Copyright (c) 2020 by Contributors
 * @file int32.H
 * @brief 32bit integer
 *
 * @author Wang Ding
 */
#ifndef __INT32_H__
#define __INT32_H__

#include <cstdint>
#include <climits>
#include <cstdlib>
#include <string>

#include "base/traits.H"

namespace Foam {

class Istream;
class Ostream;

} // namespace Foam

namespace base {

using Foam::Istream;
using Foam::Ostream;

// toString() and IOstream operators
std::string toString(const int32_t val);

int32_t readInt32(Istream& is);

bool read(const char* buffer, int32_t& val);

Istream& operator>> (Istream& is, int32_t& val);

Ostream& operator<< (Ostream& os, const int32_t val);

/*!
 * @note On 32bit OSs, long is ambiguously causing problems for IO operator resolution.
 *       This problem is avoided by explicitly defining the following operators.
 */
#if WM_ARCH_OPTION == 32
  Istream& operator>> (Istream&, long&);
  Ostream& operator<< (Ostream&, const long);
#endif

//! @brief Template specialization for pTraits<int32_t>
template<>
class Traits<int32_t> {

public:
  typedef int32_t pType;

  static const char* const TypeName;
  static const int32_t Zero;
  static const int32_t One;
  static const int32_t Min;
  static const int32_t Max;

public:

  //! @brief Construct from int32_t explicitly.
  explicit Traits(const int32_t& p): p_(p) {}

  //! @brief Construct from Istream.
  Traits(Foam::Istream& is) { is >> p_; }

  //! @brief Access to int32_t value
  operator int32_t() const { return p_; }
  operator int32_t&() { return p_; }

private:

  int32_t p_;
};

const int32_t Traits<int32_t>::Zero = 0;
const int32_t Traits<int32_t>::One = 1;
const int32_t Traits<int32_t>::Min = INT32_MIN;
const int32_t Traits<int32_t>::Max = INT32_MAX;

} // namespace base

#endif // __INT32_H__
