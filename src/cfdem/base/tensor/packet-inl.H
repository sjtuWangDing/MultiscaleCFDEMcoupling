#ifndef __PACKET_INL_H__
#define __PACKET_INL_H__

#include "logging.H"

namespace base {
namespace packet {

enum PacketArch {
  kPlain,
  kSSE2
};

#if CFDEM_USE_SSE
#define CFDEM_DEFAULT_PACKET base::packet::kSSE2
#else
#define CFDEM_DEFAULT_PACKET base::packet::kPlain
#endif

template<typename DType, PacketArch Arch = CFDEM_DEFAULT_PACKET>
struct Packet;

} // namespace packet
} // namespace base

#include "packet_plain-inl.H"
#if CFDEM_USE_SSE
#include "packet_sse-inl.H"
#endif

#endif // __PACKET_INL_H__
