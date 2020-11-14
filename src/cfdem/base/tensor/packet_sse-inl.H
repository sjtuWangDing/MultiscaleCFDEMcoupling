#ifndef __PACKET_SSE_INL_H__
#define __PACKET_SSE_INL_H__

#include <emmintrin.h> // for use sse2
#include "packet-inl.H"

namespace base {
namespace packet {

/*!
 * \brief vector real type for float
 * \note 使用 Intel C++ Compiler 或是 Microsoft Visual C++ 中提供的支持SSE指令集的 intrinsics 内联函数
 *       intrinsics 内联函数是对 MMX、SSE 等指令集的一种封装，在编译的时候，这些函数会被内联为汇编，不会产生函数调用的开销
 */
template<>
struct Packet<float, kSSE2> {
  /*! \brief number of float in vector */
  static constexpr index_t size = 4;
  /*!
   * \brief __m128 单精度浮点数，定义为：
   *        typedef struct __declspec(intrin_type) __declspec(align(16)) __m128 {
   *          float m128_f32[4]; // __declspec(align(16)) 表示该数组是以16字节为边界对齐的
   *        } __m128;
   * \note sizeof(__m128) = 16
   */
  __m128 data_;
  Packet() {}
  explicit Packet(__m128 data): data_(data) {}
  /*!
   * \brief 从内存加载数据到寄存器，将四个单精度浮点数加载到寄存器
   *        r0 := p[0], r1 := p[1], r2 := p[2], r3 := p[3]
   * \note __m128 _mm_load_ss(float *p)
   *       将一个单精度浮点数加载到寄存器的第一个字节，其它三个字节清零，即 r0 := *p, r1 := r2 := r3 := 0.0
   */
  inline static Packet<float, kSSE2> load(const float* src) {
    return Packet<float, kSSE2>(_mm_load_ps(src));
  }
  /*! \brief 从内存加载数据到寄存器，与 _mm_load_ps 一样的加载，但是不要求地址是16字节对齐 */
  inline static Packet<float, kSSE2> loadUnAligned(const float* src) {
    return Packet<float, kSSE2>(_mm_loadu_ps(src));
  }
  /*! \brief 与 load 类似，不需要字节对齐，需要多条指令 r0 = r1 = r2 = r3 = s */
  inline static Packet<float, kSSE2> fill(float s) {
    return Packet<float, kSSE2>(_mm_set1_ps(s));
  }
  inline Packet<float, kSSE2>& operator= (float s) {
    data_ = _mm_set1_ps(s);
    return *this;
  }
  /*! \brief 将 SSE 寄存器中的数据保存到内存中，dest[i] = data_[i] */
  inline void Store(float* dest) const {
    _mm_store_ps(dest, data_);
  }
};

/*! \brief vector real type for double */
template<>
struct Packet<double, kSSE2> {
  /*! \brief number of double in vector */
  static constexpr index_t size = 2;
  __m128d data_;
  Packet(void) {}
  explicit Packet(__m128d data) : data_(data) {}
  inline static Packet<double, kSSE2> load(const double* src) {
    return Packet<double, kSSE2>(_mm_load_pd(src));
  }
  inline static Packet<double, kSSE2> loadUnAligned(const double* src) {
    return Packet<double, kSSE2>(_mm_loadu_pd(src));
  }
  inline static Packet<double, kSSE2> fill(double d) {
    return Packet<double, kSSE2>(_mm_set1_pd(d));
  }
  inline Packet<double, kSSE2>& operator= (double d) {
    data_ = _mm_set1_pd(d);
    return *this;
  }
  inline void store(double *dest) const {
    _mm_store_pd(dest, data_);
  }
};

inline Packet<float, kSSE2> operator+ (const Packet<float, kSSE2>& lhs,
                                       const Packet<float, kSSE2>& rhs) {
  return Packet<float, kSSE2>(_mm_add_ps(lhs.data_, rhs.data_));
}

inline Packet<double, kSSE2> operator+ (const Packet<double, kSSE2>& lhs,
                                        const Packet<double, kSSE2>& rhs) {
  return Packet<double, kSSE2>(_mm_add_pd(lhs.data_, rhs.data_));
}

inline Packet<float, kSSE2> operator- (const Packet<float, kSSE2>& lhs,
                                       const Packet<float, kSSE2>& rhs) {
  return Packet<float, kSSE2>(_mm_sub_ps(lhs.data_, rhs.data_));
}

inline Packet<double, kSSE2> operator- (const Packet<double, kSSE2>& lhs,
                                        const Packet<double, kSSE2>& rhs) {
  return Packet<double, kSSE2>(_mm_sub_pd(lhs.data_, rhs.data_));
}

inline Packet<float, kSSE2> operator* (const Packet<float, kSSE2>& lhs,
                                       const Packet<float, kSSE2>& rhs) {
  return Packet<float, kSSE2>(_mm_mul_ps(lhs.data_, rhs.data_));
}

inline Packet<double, kSSE2> operator* (const Packet<double, kSSE2>& lhs,
                                        const Packet<double, kSSE2>& rhs) {
  return Packet<double, kSSE2>(_mm_mul_pd(lhs.data_, rhs.data_));
}

inline Packet<float, kSSE2> operator/ (const Packet<float, kSSE2>& lhs,
                                       const Packet<float, kSSE2>& rhs) {
  return Packet<float, kSSE2>(_mm_div_ps(lhs.data_, rhs.data_));
}

inline Packet<double, kSSE2> operator/ (const Packet<double, kSSE2>& lhs,
                                        const Packet<double, kSSE2>& rhs) {
  return Packet<double, kSSE2>(_mm_div_pd(lhs.data_, rhs.data_));
}

} // namespace packet
} // namespace base

#endif // __PACKET_SSE_INL_H__
