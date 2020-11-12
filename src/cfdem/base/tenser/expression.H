#ifndef __EXPRESSION_H__
#define __EXPRESSION_H__

#include <iostream>

/*!
 * \brief type of expressions defined as bitmask
 *   subtype relationshop kRValue < kMapper < kPull < kComplex
 */
namespace type {
/*!
 * \brief this expression directly correspnds to a data class,
 *   can be used to assign data
 */
const int kRValue = 0;
/*!
 * \brief expression contains element-wise tensor operations,
 *   map a expression to same shape
 */
const int kMapper = 1;

}  // namespace type

/*! \brief namespace for savers */
namespace sv {

struct SaveTo {
  /*! \brief save b to a using save method */
  template<typename DType>
  inline static void save(DType& a, DType b) { a = b; }
};

} // namespace sv

template<typename Saver, typename 

/*!
 * \brief expression engine that actually interprets these expressions
 *   this is a function template that needed to be implemented for specific expressions
 * \tparam Saver the save method
 * \tparam RValue the type of RValue to be saved
 */
template<typename Saver, typename RValue, typename DType>
struct ExpEngine {
  template<typename E>
  inline static void evaluate(RValue *dst, Exp<E, DType, type::kMapper>& exp) {
    MapExp<Saver>(dst, exp);
  }
};

/*!
 * \brief base class for expression
 * \tparam SubType inheritated class must put their type into this parameter
 * \tparam DType the data type of each element in the expression
 * \tparam exp_type expression type, see namespace type
 */
template<typename SubType, typename DType, int exp_type>
struct Exp {
 public:
  /*! \return  subtype instance of current class */
  inline const SubType& self(void) const {
    return *static_cast<const SubType*>(this);
  }
  /*! \return reference of subtype instance of current class */
  inline SubType* ptrself(void) {
    return static_cast<SubType*>(this);
  }
};

/*!
 * \brief scalar expression
 * \tparam DType the data type of the scalar
 */
template<typename DType>
struct ScalarExp: public Exp<ScalarExp<DType>, DType, type::kMapper> {
  /*! \brief scalar value */
  DType scalar_;
  /*! \brief implicit constructor, MUST NOT BE explicit */
  ScalarExp(DType scalar) : scalar_(scalar) {}  // NOLINT(*)
};

/*! \brief create an scalar expression */
template<typename DType>
inline ScalarExp<DType> scalar(DType s) {
  return ScalarExp<DType>(s);
}

/*!
 * \brief base class of all rvalues
 * \tparam Container the actually class of data container, e.g. Tensor1D
 * \tparam DataType the element data type of each element in the container
 */
template<typename Container, typename DType>
struct RValueExp: public Exp<Container, DType, type::kRValue> {
  inline Container& __assign(DType s) {
    ExpEngine<sv::SaveTo, Container, DType>::evaluate(this->ptrself(), scalar<Dtype>(s));
    return *(this->ptrself());
  }
};

template<typename Container, int dimension, typename DType>
struct TRValue: public RValueExp<> {

};

#endif // __EXPRESSION_H__
