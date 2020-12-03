#ifndef __X_ALLOCATOR_H__
#define __X_ALLOCATOR_H__

#include <new>
#include <climits> // UINT_MAX

namespace base {

/*!
 * \brief 分配内存空间（throw std::bad_alloc）
 * \tparam _T 元素类型
 * \param _N 元素个数
 * \param _p 仅仅用于函数模板的类型推导，本身在函数中没有作用
 */
template<typename _T>
inline _T* _allocate(ptrdiff_t _N, _T* _p) throw(std::bad_alloc) {
  // 卸载目前的内存分配异常处理函数，设置为系统默认函数，即 throw std::bad_alloc
  std::set_new_handler((std::new_handler)0);
  // 调用 ::operator new 分配内存空间
  _T* ptr = (_T)::operator new((size_t)(_N * sizeof(_T)));
  if (0 == ptr) {
    std::cerr << "out of memory" << std::endl;
    exit(1);
  }
  return ptr;
}

/*!
 * \brief 释放内存空间
 * \tparam _T 元素类型
 */
template<typename _T>
inline void _deallocate(_T* _buffer) throw() {
  // 必须调用 ::operator delete 释放由 ::operator new 申请的空间
  ::operator delete((void*)_buffer);
}

/*!
 * \brief base::allocator 类似于 std::allocator，可以用于替换标准容器模板中最后一个模板参数
 *   Eg: std::vector<int, base::allocator>
 *   allocator 能够提供更好的性能和更灵活的内存管理能力：
 *     - 将内存分配与对象构造分开
 *     - 将内存释放与对象析构分开
 *     - std::allocator 的实现没有成员变量，allocator 是无状态的，所有 allocator 实例都是可以互相替换的
 *     Eg:
 *     std::allocator<int> a1;
 *     std::allocator<int> a2;
 *     int* p = a1.allocate(1); // 分配一个int所需要的空间
 *     a2.dealloacte(p1); // 使用 a2 释放
 */
template<typename _T>
class allocator {
public:
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  typedef _T*       pointer;
  typedef const _T* const_pointer;
  typedef _T&       reference;
  typedef const _T& const_reference;
  typedef _T        value_type;
  /*!
   * \brief 将一个 allocator 类型转换成另外一个 allocator 类型
   *   Eg:
   *   base::allocator<int> a1;
   *   decltype(a1)::rebind<char>::other a2;
   */
  template<class _T1>
  struct rebind {
    typedef allocator<_T1> other;
  };
  //! \brief Constructor
  allocator() throw() {}
  allocator(const allocator& alloc) throw() {}
  template <class _T1>
  allocator(const allocator<_T1>& alloc) throw() {}
  //! \brief Deconstructor
  ~allocator() {}
  //! \brief 取地址
  pointer address(reference _x ) const { return &_x; }
  //! \brief 取地址
  const_pointer address(const_reference _x) const { return &_x; }
  //! \brief 最大容量
  size_type max_size() throw() { return (size_type)(UINT_MAX / sizeof(_T)); }
  //! \brief 分配 _N 个 _T 类型的内存
  pointer allocate(size_type _N) {
    // 卸载目前的内存分配异常处理函数，设置为系统默认函数，即 throw std::bad_alloc
    std::set_new_handler((std::new_handler)0);
    if (_N > max_size()) {
      throw std::bad_alloc();
    }
    // 调用 ::operator new 分配内存空间
    return (pointer)::operator new((size_t)(_N * sizeof(_T)));
  }
  //! \brief 释放内存
  void deallocate(pointer _buffer) {
    // 必须调用 ::operator delete 释放由 ::operator new 申请的空间
    ::operator delete((void*)_buffer);
  }
  //! \brief 构造对象
  template<typename _U, typename... _Args>
  void construct(_U* _p, _Args&&... _args) {
    ::new ((void*)_p) _U(std::forward<_Args>(_args)...);
  }
  //! \brief 析构对象
  template<typename _U>
  void destroy(_U* _p) {
    _p->~U();
  }
};

} // namepsace base

#endif // __X_ALLOCATOR_H__
