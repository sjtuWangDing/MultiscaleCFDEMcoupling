#ifndef __TENSOR_CPU_INL_H__
#define __TENSOR_CPU_INL_H__

#include "tensor.H"

namespace base {

// template<typename Saver, typename TContainer, int dim,
//          typename DType, typename EContainer, int exp_type>
// inline void MapExp(TRValue<TContainer, dim, DType> *dst,
//                    const expr::Exp<EContainer, DType, exp_type>& exp) {
//   Shape<dim> expShape = expr::ShapeCheck<dim, EContainer>::Check(exp.self());
//   Shape<dim> destShape = expr::ShapeCheck<dim, TContainer>::Check(dst->self());
//   CHECK(expShape[0] == 0 || expShape == destShape)
//       << "Assignment: Shape of Tensors are not consistent with target, "
//       << "expression shape: " << expShape
//       << " dest shape: " << destShape;
// }

} // namespace base

#endif // __TENSOR_CPU_INL_H__
