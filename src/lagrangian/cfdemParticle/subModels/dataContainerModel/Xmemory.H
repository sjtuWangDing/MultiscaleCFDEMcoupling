/*!
 *  Copyright (c) 2020 by Contributors
 * \file Xmemory.H
 * \brief memory allocation like STL
 * \author Wang Ding
 */
#ifndef __XMEMORY_H__
#define __XMEMORY_H__

#include <new>      // placement new
#include <cstddef>  // size_t / ptrdiff_t
#include <cstdlib>  // exit()
#include <climits>  // UINT_MAX

namespace Foam {

/*!
 * \brief 分配内存空间
 * \tparam _Ty - 仅仅用于函数模板的类型推导
 */
template<typename _Ty>
inline _Ty* _Allocate(ptrdiff_t _N, _Ty*) {
	// 设置 new 处理函数为系统默认
	// 内存分配失败抛出 std::bad_alloc 异常
	std::set_new_handler((std::new_handler)NULL);
	if (_N < 0) {
		_N = 0;
	}
	// 调用 operator new 分配内存
	return (_Ty*)::operator new((size_t)_N * sizeof(_Ty));
}

/*!
 * \brief 释放内存空间
 * \tparam _Ty - 数据类型
 */
template<typename _Ty>
inline void _Deallocate(_Ty* _Buffer) {
	// 必须调用 ::operator delete 释放由 ::operator new 申请的内存
	::operator delete(_Buffer);
}

/*!
 * \brief 构造
 * \tparam _T1 - 构造对象类型
 * \tparam _T2 - 构造参数类型
 */
template<typename _T1, typename _T2>
inline void _Construct(_T1* _P, const _T2& _V) {
	// 调用 placement new 在 _P 指向的内存上构造对象
	// 这里就会调用 _T1 类中的构造函数或拷贝构造函数
	new((void*)_P) _T1(_V);
}

/*!
 * \brief 析构
 * \tparam _Ty - 对象类型
 */
template<typename _Ty>
inline void _Destroy(_Ty* _T) {
	_T->~_Ty();
}

/*!
 * \brief 空间配置器
 *        分配容器节点的内存(_Node)
 *        分配容器中每个节点中的元素的内存(_Ty)
 *        调用 _Ty 构造函数
 * \tparam _Ty - 数据类型
 */
template<typename _Ty>
class allocator {
public:
	typedef _Ty        value_type;
	typedef _Ty*       pointer;
	typedef const _Ty* const_pointer;
	typedef _Ty&       reference;
	typedef const _Ty& const_reference;
	typedef size_t     size_type;
	typedef ptrdiff_t  difference_type;
	/*!
	 * \brief allocator 重新绑定类型为 _Tu
	 * \tparam _Tu - 数据类型
	 */
	template<typename _Tu>
	struct rebind {
		typedef allocator<_Tu> other;
	};
	/*!
	 * \brief 分配节点内存
	 * \note _Node* _S = charalloc(sizeof(_Node))
	 *        分配大小为 sizeof(_Node) 的内存
	 * \param _N - 内存字节数 
	 */
	char* charalloc(size_type _N) {
		return _Allocate((difference_type)_N, (char*)NULL);
	}
	/*!
	 * \brief 分配元素内存
	 * \param _N - 元素个数
	 */
	pointer allocate(size_type _N) {
		return _Allocate((difference_type)_N, (pointer)NULL);
	}
	/*!
	 * \brief 释放内存
	 */
	void deallocate(void* _P) {
		return _Deallocate(_P);
	}
	/*!
	 * \brief 构造元素
	 * \param _P - 构造元素地址
	 * \param _V - 构造参数 
	 */
	void construct(pointer _P, const _Ty& _V) {
		_Construct(_P, _V);
	}
	/*!
	 * \brief 析构对象
	 * \param _P - 析构元素地址
	 */
	void destroy(pointer _P) {
		_Destroy(_P);
	}
	/*!
	 * \brief 取地址
	 * \param _X - 元素引用
	 */
	pointer address(reference _X) {
		return &_X;
	}
	/*!
	 * \brief 取地址
	 * \param _X - 元素常引用
	 */
	const_pointer address(const_reference _X) {
		return (const_pointer)&_X;
	}
	/*!
	 * \brief 计算容量 
	 */
	size_type max_size() const {
		return (size_type)(UINT_MAX / sizeof(value_type));
	}
};

/*!
 * \brief allocator 模板的全类特化
 */
template<>
class allocator<void> {
public:
	typedef void        value_type;
	typedef void*       pointer;
	typedef const void* const_pointer;
};

}  // End of namespace Foam

#endif  // __XMEMORY_H__
