#ifndef __TENSOR_H__
#define __TENSOR_H__

#include "logging.H"
#include "expression.H"

namespace base {

/*!
 * \brief shape of a tensor
 * \tparam dimension dimension of tensor
 */
template<int dimension>
class Shape {
 public:
  /*! \brief dimension of current shape */
  static const int kDimension = dimension;
  /*! \brief dimension of current shape minus one */
  static const int kSubdim = dimension - 1;
  /*! \brief storing the dimension information */
  index_t shape_[kDimension];
  /*! \brief default constructor, do nothing */
  inline Shape(void) {}
  /*! \brief constuctor */
  inline Shape(const Shape<kDimension>& s) {
    #pragma unroll
    for (int i = 0; i < kDimension; ++i) {
      this->shape_[i] = s.shape_[i];
    }
  }
  /*! \return number of valid elements */
  inline index_t size() const {
    index_t res = this->shape_[0];
    for (int i = 1; i < kDimension; ++i) {
      res *= this->shape_[i];
    }
    return res;
  }
  /*!
   * \brief get corresponding index
   * \param idx dimension index
   * \return the corresponding dimension size
   */
  inline index_t& operator[] (int idx) {
    CHECK(idx >= 0 && idx < kDimension);
    return shape_[idx];
  }
  inline index_t operator[] (int idx) const {
    CHECK(idx >= 0 && idx < kDimension);
    return shape_[idx];
  }
  /*!
   * \return whether two shape equals
   * \param s the shape to compare against
   */
  inline bool operator== (const Shape<kDimension> &s) const {
    #pragma unroll
    for (int i = 0; i < kDimension; ++i) {
      if (s.shape_[i] != this->shape_[i]) { return false; }
    }
    return true;
  }
  /*!
   * \return whether two shape not equal
   * \param s the shape to compare against
   */
  inline bool operator!= (const Shape<kDimension>& s) const {
    return !(*this == s);
  }
  /*!
   * \brief get subshape that takes off largest dimension
   * \return subshape
   */
  inline Shape<kSubdim> subShape() const {
    Shape<kSubdim> s;
    #pragma unroll
    for (int i = 0; i < kSubdim; ++i) {
      s.shape_[i] = this->operator[](i + 1);
    }
    return s;
  }
  /*!
   * \brief flatten the tensor, return a 1D shape
   * \return the flat 1d shape
   */
  Shape<1> flatTo1D() const {
    Shape<1> s;
    s[0] = this->size();
    return s;
  }
  /*! \brief allow string printing of the shape */
  template<int dim>
  friend std::ostream& operator<< (std::ostream& os, const Shape<dim>& shape);
};

template<int dim>
std::ostream& operator<< (std::ostream& os, const Shape<dim>& shape) {
  os << "(";
  for (int i = 0; i < dim; ++i) {
    if (i != 0) os << ", ";
    os << shape[i];
  }
  os << ")";
  return os;
}

/*!
 * \brief construct a one dimension shape, stride will equal s0
 * \param s0 size of dimension 0
 * \return the shape construction
 */
inline Shape<1> makeShape1(index_t s0) {
  Shape<1> s;
  s[0] = s0;
  return s;
}

/*!
 * \brief construct a two dimension shape, stride will equal s0
 * \param s0 size of dimension 0
 * \param s1 size of dimension 1
 * \return the shape construction
 */
inline Shape<2> makeShape2(index_t s0, index_t s1) {
  Shape<2> s;
  s[0] = s0;
  s[1] = s1;
  return s;
}

/*!
 * \brief Tensor RValue, this is the super type of all kinds of possible tensors
 * \tparam Container the tensor type
 * \tparam dimension dimension of the tensor
 * \tparam Device which device the tensor is on
 * \tparam DType the type of elements in the tensor
 */
template<typename Container, int dimension, typename Device, typename DType>
struct TRValue: public expr::RValueExp<Container, DType> {};

/*!
 * \brief general tensor
 * \tparam dimension dimension of the tensor
 * \tparam Device which device the tensor is on
 * \tparam DType the type of elements in the tensor
 */
template<int dimension, typename Device CFDEM_DEFAULT_DEVICE, typename DType CFDEM_DEFAULT_DTYPE>
class Tensor: public TRValue<Tensor<dimension, Device, DType>,
                             dimension, Device, DType> {
 public:
  /*! \brief whether current type lies in cpu */
  static const bool kDevCPU = Device::kDevCPU;
  /*! \brief dimension of subtype */
  static const int  kSubdim = dimension - 1;
  /*! \brief pointer to the data */
  DType *dptr_;
  /*! \brief shape of the tensor */
  Shape<dimension> shape_;
  /*!
   * \brief storing the stride information in x dimension
   *        this is used to deal with pitch allocation in gpu or sse(align x dimension to 64bit) for efficiency
   */
  index_t stride_;
  /*! \brief default constructor */
  inline Tensor(): dptr_(nullptr), shape_(), stride_(0) {}
  /*! \brief constructor from shape */
  inline Tensor(const Shape<dimension>& shape)
      : dptr_(nullptr), shape_(shape), stride_(shape[kSubdim]) {}
  /*! \brief constructor from data pointer and shape, without stride */
  inline Tensor(DType *dptr, const Shape<dimension>& shape)
      : dptr_(dptr), shape_(shape), stride_(shape[kSubdim]) {}
  /*! \brief constructor from data pointer and shape */
  inline Tensor(DType *dptr, const Shape<dimension>& shape, index_t stride)
      : dptr_(dptr), shape_(shape), stride_(stride) {}
  /*! \brief implement the assignment of same type */
  inline Tensor<dimension, DType>& operator= (const Tensor<dimension, DType>& exp) {
    dptr_ = exp.dptr_;
    shape_ = exp.shape_;
    stride_ = exp.stride_;
    return *this;
  }
  /*!
   * \return memory cost of the tensor, including the aligned x dimension
   * \tparam startdim the starting dimension
   */
  template<int startDim>
  inline index_t memSize() const {
    if (dptr_ == nullptr) { return 0; }
    index_t res = this->stride_;
    #pragma unroll
    for (int i = startDim; i < kSubdim; ++i) {
      res *= this->shape_[i];
    }
    return res;
  }
  /*! \return whether the tensor's memory is continuous */
  inline bool isContinus() const {
    return this->shape_[kSubdim] == stride_;
  }
  /*! \return memory cost of the tensor, including the aligned x dimension */
  inline index_t mSize() const {
    return this->memSize<0>();
  }
  /*!
   * \brief return size of i-th dimension, start counting from highest dimension
   * \param idx the dimension count from the highest dimensin
   * \return the size
   */
  inline index_t size(int idx) const {
    return shape_[idx];
  }
  /*!
   * \brief get a element of dimension - 1
   * \param idx index
   * \return the result tensor
   */
  inline Tensor<kSubdim, Device, DType> operator[] (int idx) const {
    return Tensor<kSubdim, Device, DType>(dptr_ + this->memSize<1>() * idx, shape_.subShape(), stride_);
  }
};

/*! \brief respecialized class Tensor1D, thei is due to different implementation in operator[] */
template<typename Device, typename DType>
class Tensor<1, Device, DType>: public TRValue<Tensor<1, Device, DType>,
                                               1, Device, DType> {
 public:
  DType *dptr_;
  Shape<1> shape_;
  index_t stride_;
  // constructors
  inline Tensor(): dptr_(nullptr), shape_(), stride_(0) {}
  inline Tensor(const Shape<1>& shape)
      : dptr_(nullptr), shape_(shape), stride_(shape[0]) {}
  inline Tensor(DType *dptr, const Shape<1> &shape)
      : dptr_(dptr), shape_(shape), stride_(shape[0]) {}
  inline Tensor(DType *dptr, const Shape<1> &shape, index_t stride)
      : dptr_(dptr), shape_(shape), stride_(stride) {}
  inline Tensor<1, Device, DType>& operator= (const Tensor<1, Device, DType> &exp) {
    dptr_ = exp.dptr_;
    shape_ = exp.shape_;
    stride_ = exp.stride_;
    return *this;
  }
  // member function
  template<int startDim>
  inline index_t memSize() const {
    CHECK_EQ(startDim, 0);
    return stride_;
  }
  inline bool isContinus() const {
    return stride_ == shape_[0];
  }
  inline index_t mSize() const {
    return shape_[0];
  }
  inline index_t size(int idx) const {
    return shape_[idx];
  }
  inline DType& operator[] (int idx) const {
    CHECK(idx >= 0 && idx < static_cast<int>(mSize()));
    return dptr_[idx];
  }
};

using CITensor1 = ::base::Tensor<1, cpu, int>;
using CDTensor1 = ::base::Tensor<1, cpu, double>;
using CITensor2 = ::base::Tensor<2, cpu, int>;
using CDTensor2 = ::base::Tensor<2>;

//------------------------
// Function Declarations
//-----------------------

// /*!
//  * \brief map a expression to a tensor, this function calls MapPlan
//  *        implemented in tensor_cpu-inl.H
//  * \tparam Saver specify storage method
//  * \tparam TContainer specifies the storage type of the tensor
//  * \tparam dim dim of the tensor, during usage, there is no need to specify this parameter
//  * \tparam DType the type of elements in the tensor
//  * \tparam EContainer specifies the expression type, not need to specify this parameter during usage
//  * \tparam exp_type expression type
//  * \param dst destination
//  * \param exp expression
//  */
// template<typename Saver, typename TContainer, int dim,
//          typename DType, typename EContainer, int exp_type>
// inline void MapExp(TRValue<TContainer, dim, DType> *dst,
//                    const expr::Exp<EContainer, DType, exp_type>& exp);
// {
//   Shape<dim> expShape = expr::ShapeCheck<dim, EContainer>::Check(exp.self());
//   Shape<dim> destShape = expr::ShapeCheck<dim, TContainer>::Check(dst->self());
//   CHECK(expShape[0] == 0 || expShape == destShape)
//       << "Assignment: Shape of Tensors are not consistent with target, "
//       << "expression shape: " << expShape
//       << " dest shape: " << destShape;
// }

} // namespace base

#include "expression_engine-inl.H"
#include "tensor_cpu-inl.H"

#endif // __TENSOR_H__
