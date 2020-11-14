#ifndef __PACKET_PLAIN_INL_H__
#define __PACKET_PLAIN_INL_H__

#include "packet-inl.H"

namespace base {
namespace packet {

template<typename DType>
struct Packet<DType, kPlain> {
 public:
  /*! \brief number of float in vector */
  static constexpr index_t size = 1;
  /*! \brief The internal data */
  DType data_;
  /*! \brief Default constructor */
  Packet(void) {}
  /*! \brief Construct from data */
  explicit Packet(DType data): data_(data) {}
  /*! \brief Overload operator= */
  Packet<DType, kPlain>& operator= (DType s) {
    data_ = s;
    return *this;
  }
  /*! \brief Overload operator= */
  Packet<DType, kPlain>& operator= (const Packet<DType, kPlain>& ps) {
    data_ = ps.data_;
    return *this;
  }
  /*! \brief Fill with s */
  inline static Packet<DType, kPlain> fill(DType s) {
    return Packet<DType, kPlain>(s);
  }
  /*! \brief load from address */
  inline static Packet<DType, kPlain> load(const DType* src) {
    return Packet<DType, kPlain>(*src);
  }
  /*! \brief load from un aligned address */
  inline static Packet<DType, kPlain> loadUnAligned(const DType* src) {
    return Packet<DType, kPlain>(*src);
  }
  /*! \brief store data to address */
  inline void store(DType* dest) const {
    *dest = data_;
  }
};

template<typename DType>
inline Packet<DType, kPlain> operator+ (const Packet<DType, kPlain>& lhs,
                                        const Packet<DType, kPlain>& rhs) {
  return Packet<DType, kPlain>(lhs->data_ + rhs->data_);
}

template<typename DType>
inline Packet<DType, kPlain> operator- (const Packet<DType, kPlain>& lhs,
                                        const Packet<DType, kPlain>& rhs) {
  return Packet<DType, kPlain>(lhs->data_ - rhs->data_);
}

template<typename DType>
inline Packet<DType, kPlain> operator* (const Packet<DType, kPlain>& lhs,
                                        const Packet<DType, kPlain>& rhs) {
  return Packet<DType, kPlain>(lhs->data_ * rhs->data_);
}

template<typename DType>
inline Packet<DType, kPlain> operator/ (const Packet<DType, kPlain>& lhs,
                                        const Packet<DType, kPlain>& rhs) {
  return Packet<DType, kPlain>(lhs->data_ / rhs->data_);
}

} // namespace packet
} // namespace base

#endif // __PACKET_PLAIN_INL_H__
