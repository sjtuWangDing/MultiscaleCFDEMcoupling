/*---------------------------------------------------------------------------*\
  CFDEMcoupling - Open Source CFD-DEM coupling

  CFDEMcoupling is part of the CFDEMproject
  www.cfdem.com
                              Christoph Goniva, christoph.goniva@cfdem.com
                              Copyright 2009-2012 JKU Linz
                              Copyright 2012-     DCS Computing GmbH, Linz
------------------------------------------------------------------------------
License
  This file is part of CFDEMcoupling.

  CFDEMcoupling is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 3 of the License, or (at your
  option) any later version.

  CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with CFDEMcoupling; if not, write to the Free Software Foundation,
  Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
\*---------------------------------------------------------------------------*/

#ifndef __RUN_TIME_SELECTION_TABLES_H__
#define __RUN_TIME_SELECTION_TABLES_H__

#include <unordered_map>
#include "error.H"
#include "word.H"

/*!
 * @brief 在所有 sub model 的类声明中声明 type name 以及函数
 *        无论是这个模型的基类，还是子类都需要声明
 * @param typeNameCString type name，类型应该是 const char*
 */
#define cfdemTypeName(typeNameCString)                                                                             \
  static const char* cTypeName() { return typeNameCString; }                                                       \
  static const std::string typeName_;                                                                              \
  virtual const std::string& typeName() const { return typeName_; }

/*!
 * @brief 在所有 sub model 的类定义中初始化静态成员 typeName_
 *        无论是这个模型的基类，还是子类都需要声明
 * @param Type 指定当前类的类名
 */
#define cfdemDefineTypeName(Type)                                                                                  \
  const std::string Type::typeName_(Type::cTypeName());

/*!
 * @brief 运行时选择器，这个宏只需要在每一个 sub model 的基类的声明中使用
 * @param AutoPtr 用于封装 sub model 指针的模板，通常就选择 Foam::autoPtr
 * @param BaseType 基类类型
 * @param argList 构造函数的形参
 * @param parList 构造函数的实参
 * @pre NewFunctionPtr: 定义一个函数指针类型，返回值类型为 autoPtr<baseType>，参数为 argList
 * @pre NewFunctionMap: 定义 unordered_map 的类型，key 类型为 std::string，value 类型为 NewFunctionPtr
 * @pre newFunctionMapPtr_: 声明 NewFunctionMap 类型静态指针
 * @pre constructNewFunctionMap: 创建 newFunctionMapPtr_ 的静态方法
 * @pre destroyNewFunctionMap: 析构 newFunctionMapPtr_ 的静态方法
 * @pre AddNewFunctionToMap: 在子类的实现中创建 AddNewFunctionToMap<子类> 的一个实例
 *                           则这个实例在构造方法中，会将 AddNewFunctionToMap<子类>::New 这个静态方法插入到 map 中，
 *                           而这个静态方法就是实际用来创建子类对象的方法
 */
#define cfdemDeclareRunTimeSelection(AutoPtr, BaseType, argList, parList)                                          \
                                                                                                                   \
    typedef AutoPtr<BaseType> (*NewFunctionPtr)argList;                                                            \
                                                                                                                   \
    typedef std::unordered_map<std::string, NewFunctionPtr> NewFunctionMap;                                        \
                                                                                                                   \
    static NewFunctionMap* newFunctionMapPtr_;                                                                     \
                                                                                                                   \
    static void constructNewFunctionMap();                                                                         \
                                                                                                                   \
    static void destroyNewFunctionMap();                                                                           \
                                                                                                                   \
    template<class BaseType##Type>                                                                                 \
    class AddNewFunctionToMap {                                                                                    \
    public:                                                                                                        \
      static AutoPtr<BaseType> New argList {                                                                       \
        return AutoPtr<BaseType>(new BaseType##Type parList);                                                      \
      }                                                                                                            \
      AddNewFunctionToMap(const std::string& typeName = BaseType##Type::typeName_) {                               \
        constructNewFunctionMap();                                                                                 \
        if (newFunctionMapPtr_->find(typeName) == newFunctionMapPtr_->end()) {                                     \
          newFunctionMapPtr_->emplace(typeName, New);                                                              \
        } else {                                                                                                   \
          std::cerr << "Duplicate entry: "                                                                         \
            << typeName << " in runtime selection table "                                                          \
            << #BaseType << std::endl;                                                                             \
          error::safePrintStack(std::cerr);                                                                        \
        }                                                                                                          \
      }                                                                                                            \
      ~AddNewFunctionToMap() {                                                                                     \
        destroyNewFunctionMap();                                                                                   \
      }                                                                                                            \
    };

/*!
 * @brief 定义 BaseType 中的静态方法 constructNewFunctionMap()，用于创建对象 map 对象
 * @param BaseType 指定基类
 */
#define cfdemDefineConstructNewFunctionMap(BaseType)                                                               \
  void BaseType::constructNewFunctionMap() {                                                                       \
    static bool constructed = false;                                                                               \
    if (!constructed) {                                                                                            \
      constructed = true;                                                                                          \
      BaseType::newFunctionMapPtr_ = new BaseType::NewFunctionMap;                                                 \
    }                                                                                                              \
  }

/*!
 * @brief 定义 BaseType 中的静态方法 newRunTimeSelectionConstructorMap()，用于创建对象 map 对象
 * @param BaseType 指定基类
 */
#define cfdemDefineDestroyNewFunctionMap(BaseType)                                                                 \
  void BaseType::destroyNewFunctionMap() {                                                                         \
    if (BaseType::newFunctionMapPtr_ != nullptr) {                                                                 \
      delete BaseType::newFunctionMapPtr_;                                                                         \
      BaseType::newFunctionMapPtr_ = nullptr;                                                                      \
    }                                                                                                              \
  }

/*!
 * @brief 用于在基类中定义 New 这个静态方法，而这个方法返回的就是运行期间动态选择的某一个子类对象的autoPtr
 * @param AutoPtr 用于封装 sub model 指针的模板，通常就选择 Foam::autoPtr
 * @param BaseType 指定基类
 * @param argList 构造函数的形参
 * @param parList 构造函数的实参
 * @param parCloud parList 中 cfdemCloud& 类型的实参
 */
#define cfdmeDefineBaseTypeNew(AutoPtr, BaseType, argList, parList, parCloud)                                      \
  AutoPtr<BaseType> BaseType::New argList {                                                                        \
    const dictionary& dict = parCloud.couplingPropertiesDict();                                                    \
    std::string modelType = dict.found(#BaseType) ? word(dict.lookup(#BaseType)) : "";                             \
    Info << "Selecting " << #BaseType << ": " << modelType << endl;                                                \
    if (BaseType::newFunctionMapPtr_->find(modelType) == BaseType::newFunctionMapPtr_->end()) {                    \
      FatalError << #BaseType << "::New" << #argList << ": "                                                       \
        << "unknow " << #BaseType << " type " << modelType << ", constructor not in map" << endl                   \
        << "Valid " << #BaseType << " types are: [";                                                               \
      for (const std::pair<std::string, BaseType::NewFunctionPtr>& item : *(BaseType::newFunctionMapPtr_)) {       \
        Info << item.first << ", ";                                                                                \
      }                                                                                                            \
      Info << "]" << endl << abort(FatalError);                                                                    \
    }                                                                                                              \
    BaseType::NewFunctionPtr& selectNewFunction = (BaseType::newFunctionMapPtr_)->operator[](modelType);           \
    return AutoPtr<BaseType>(selectNewFunction parList);                                                           \
  }

#define cfdemAddToNewFunctionMap(BaseType, DerivedType)                                                            \
  BaseType::AddNewFunctionToMap<DerivedType> add##DerivedType##NewFunctionToMap_();

#endif // __RUN_TIME_SELECTION_TABLES_H__
