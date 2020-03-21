/*!
 *  Copyright (c) 2020 by Contributors
 * \file Iterator.H
 * \brief Iterator like STL
 * \author Wang Ding
 */

#ifndef __ITERATOR_H__
#define __ITERATOR_H__

#include <cstddef>  // ptrdiff_t

namespace Foam {

/*!
 * \brief STL中规定的五种迭代器，其遵循概念强化(concept refinement)关系
 */
struct input_iterator_tag {};                                             // 只读迭代器
struct output_iterator_tag {};                                            // 只写迭代器
struct forward_iterator_tag: public input_iterator_tag {};                // 读写迭代器
struct bidirectional_iterator_tag: public forward_iterator_tag {};        // 双向迭代器
struct random_access_iterator_tag: public bidirectional_iterator_tag {};  // 随机迭代器

/*!
 * \brief STL规定的标准迭代器
 * \note 自定义迭代器最好继承标准迭代器以保证自定义迭代器一定包含上面五种内嵌类型
 *       类模板的类型参数可以带有缺省值(遵循向右对齐即可)，同事缺省值可以取自前面的类型形参
 */
template<typename Category, typename T, typename Distance = ptrdiff_t,
	     typename Pointer = T*, typename Reference = T&>
struct Iterator {
	typedef Category  iterator_category;  // 迭代器类型 
	typedef T         value_type;         // 迭代器所指对象类型
	typedef Distance  difference_type;    // 迭代器间距类型(or容器容量类型)
	typedef Pointer   pointer;            // 迭代器所指对象的指针类型
	typedef Reference reference;          // 迭代器所指对象的引用类型
};

/*!
 * \brief 迭代器类型萃取器
 * \note 需要使用 typename 解决嵌套依赖
 */
template<typename Iterator>
struct iterator_traits {
	typedef typename Iterator::iterator_category iterator_category;
	typedef typename Iterator::value_type        value_type;
	typedef typename Iterator::difference_type   difference_type;
	typedef typename Iterator::pointer           pointer;
	typedef typename Iterator::reference         reference;
};

/*!
 * \brief 迭代器类型萃取器
 * \note 针对指针(native pointer)模板偏特化 
 */
template<typename T>
struct iterator_traits<T*> {
	typedef random_access_iterator_tag iterator_category;  // 指针本质上就是随机访问迭代器
	typedef T                          value_type;
	typedef ptrdiff_t                  difference_type;
	typedef T*                         pointer;
	typedef T&                         reference;
};

/*!
 * \brief 迭代器类型萃取器
 * \note 针对常量指针(native const pointer)模板偏特化 
 */
template<typename T>
struct iterator_traits<const T*> {
	typedef random_access_iterator_tag iterator_category;  // 指针本质上就是随机访问迭代器
	typedef T                          value_type;
	typedef ptrdiff_t                  difference_type;
	typedef T*                         pointer;
	typedef T&                         reference;
};

/*!
 * \brief 获取迭代器自身类型 
 */
template<typename Iterator>
inline typename iterator_traits<Iterator>::iterator_category
iterator_category(const Iterator&) {
	typedef typename iterator_traits<Iterator>::iterator_category category;
	return category();
}

/*! \brief 获取迭代器距离类型 */
template<typename Iterator>
inline typename iterator_traits<Iterator>::difference_type*
distance_type(const Iterator&) {
	return static_cast<typename iterator_traits<Iterator>::difference_type*>(NULL);
}

/*!
 * \brief 获取迭代器指向元素的类型
 */
template<typename Iterator>
inline typename iterator_traits<Iterator>::value_type*
value_type(const Iterator&) {
	return static_cast<typename iterator_traits<Iterator>::value_type*>(NULL);
}

/*!
 * \brief 计算两个迭代器之间的距离(内部调用)
 * \tparam InputIterator - 类型参数是 input_iterator_tag 类或者其子类
 * \param first - 起始迭代器
 * \param last - 终止迭代器
 * \param input_iterator_tag - 迭代器类型(通过函数参数实现重载)
 */
template<typename InputIterator>
inline typename iterator_traits<InputIterator>::difference_type
__distance(InputIterator first,
					 InputIterator last,
					 input_iterator_tag) {
	typename iterator_traits<InputIterator>::difference_type n = 0;
	while (first != last) {
		++first;
		++n;
	}
	return n;
}

/*!
 * \brief 计算两个迭代器之间的距离(内部调用)
 * \tparam InputIterator - 类型参数是 random_access_iterator_tag
 * \param first - 起始迭代器
 * \param last - 终止迭代器
 * \param input_iterator_tag - 迭代器类型(通过函数参数实现重载)
 */
template<typename RandomAccessIterator>
inline typename iterator_traits<RandomAccessIterator>::difference_type
__distance(RandomAccessIterator first,
					 RandomAccessIterator last,
					 RandomAccessIterator) {
	return last - first;
}

/*! \brief 计算两个迭代器之间的距离(外部调用) */
template<typename InputIterator>
inline typename iterator_traits<InputIterator>::difference_type
distance(InputIterator first, InputIterator last) {
	// 获取迭代器类型
	typedef typename iterator_traits<InputIterator>::iterator_category category;
	// 利用迭代器类型构造临时对象，实现 __distance 函数的重载选择
	return __distance(first, last, category());
}

/*!
 * \brief 移动迭代器(内部调用)
 *        将迭代器 i 移动 +n 位，迭代器类型参数是 input_iterator_tag 以及其子类
 */
template<typename InputIterator, typename Distance>
inline void __advance(InputIterator& i, Distance n, input_iterator_tag) {
	while (n--) {
		++i;
	}
}

/*!
 * \brief 移动迭代器(内部调用)
 *        将迭代器 i 移动 +n or -n 位，迭代器类型参数是 bidirectional_iterator_tag 以及其子类
 */
template<typename BidirectionalIterator, typename Distance>
inline void __advance(BidirectionalIterator& i,
											Distance n,
											bidirectional_iterator_tag) {
	if (n >= 0) {
		while (n--) {
			++i;
		}
	} else {
		while (n++) {
			--i;
		}
	}
}

/*!
 * \brief 移动迭代器(内部调用)
 *        将迭代器 i 移动 +n or -n 位，迭代器类型参数是 random_access_iterator_tag 以及其子类
 */
template<typename RandomAccessIterator, typename Distance>
inline void __advance(RandomAccessIterator& i,
											Distance n,
											random_access_iterator_tag) {
	i += n;
}

/*! \brief 移动迭代器 */
template<typename InputIterator, typename Distance>
inline void advance(InputIterator& i, Distance n) {
	__advance(i, n, iterator_category(i));
}

}  // End of namespace Foam

#endif  // __ITERATOR_H__
