#ifndef __EXPRESSION_ENGINE_INL_H__
#define __EXPRESSION_ENGINE_INL_H__

#include "expression.H"
#include "tensor.H"

namespace base {
namespace expr {

/*!
 * \brief runtime shape checking template
 *    get the shape of an expression, report error if shape mismatch
 * \tparam dim the dimension of the shape
 * \tparam E expression
 */
template<int dim, typename E>
struct ShapeCheck {
  inline static Shape<dim> check(const E& t);
};

/*! \brief specialized class for ScalarExp */
template<int dim, typename DType>
struct ShapeCheck<dim, ScalarExp<DType> > {
  inline static Shape<dim> check(const ScalarExp<DType>& exp) {
    Shape<dim> shape;
    for (int i = 0; i < dim; ++i) {
      shape[i] = 0;
    }
    return shape;
  }
};

/*! \brief specialized class for Tensor */
template<int dim, typename DType>
struct ShapeCheck<dim, Tensor<dim, DType> > {
  inline static Shape<dim> check(const Tensor<dim, DType>& exp) {
    return exp.shape_;
  }
};

// /*!
//  * \brief expression engine that actually interprets these expressions
//  *   this is a function template that needed to be implemented for specific expressions
//  * \tparam Saver the save method
//  * \tparam RValue the type of RValue to be saved
//  */
// template<typename Saver, typename RValue, typename DType>
// struct ExpEngine {
//   template<typename E>
//   inline static void evaluate(RValue *dst, Exp<E, DType, type::kRValue>& exp) {
//     MapExp<Saver>(dst, exp);
//   }
//   template<typename E>
//   inline static void evaluate(RValue *dst, Exp<E, DType, type::kMapper>& exp) {
//     MapExp<Saver>(dst, exp);
//   }
// };

} // namespace expr
} // namespace base

#endif // __EXPRESSION_ENGINE_INL_H__
