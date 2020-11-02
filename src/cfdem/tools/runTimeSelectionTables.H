#ifndef __RUN_TIME_SELECTION_TABLES_H__
#define __RUN_TIME_SELECTION_TABLES_H__

#include "autoPtr.H"
#include "HashTable.H"

namespace cfdem {

// template<class T, class Key, class Hash> class HashTable;
#define declareRunTimeSelectionTable(autoPtr, baseType, argNames, argList, parList) \
  typedef autoPtr<baseType> (*argNames##ConstructorPtr)argList; \
  typedef HashTable<argNames##ConstructorPtr, word, string::hash> argNames##ConstructorTable; \
  static argNames##ConstructorTable* argNames##ConstructorTablePtr_; \
  static void construct##argNames##ConstructorTables(); \
  static void destroy##argNames##ConstructorTables(); \
  \
  template<class baseType##Type> \
  class add##argNames##ConstructorToTable \
  { \
  public: \
    static autoPtr<baseType> New argList \
    { \
      return autoPtr<baseType>(new baseType##Type parList); \
    } \
    \
    add##argNames##ConstructorToTable \
    ( \
      const word& lookup = baseType##Type::typeName \
    ) \
    { \
      construct##argNames##ConstructorTables(); \
      if (!argNames##ConstructorTablePtr_->insert(lookup, New))          \
      {                                                                  \
        std::cerr << "Duplicate entry " << lookup                       \
          << " in runtime selection table " << #baseType             \
          << std::endl;                                              \
        error::safePrintStack(std::cerr);                              \
      }                                                                  \
    } \
    ~add##argNames##ConstructorToTable()                                   \
    {                                                                      \
      destroy##argNames##ConstructorTables();                            \
    }                                                                      \
  }

#define defineRunTimeSelectionTableConstructor(baseType, argNames) \
  void baseType::construct##argNames##ConstructorTables() { \
    static bool constructed = false; \
    if (!constructed) { \
      constructed = true; \
      baseType::argNames##ConstructorTablePtr_ = new baseType::argNames##ConstructorTable; \
    } \
  }

#define defineRunTimeSelectionTableDestructor(baseType,argNames)               \
                                                                               \
    /* Table destructor called from the table add function destructor */       \
    void baseType::destroy##argNames##ConstructorTables()                      \
    {                                                                          \
      if (baseType::argNames##ConstructorTablePtr_)                          \
      {                                                                      \
        delete baseType::argNames##ConstructorTablePtr_;                   \
        baseType::argNames##ConstructorTablePtr_ = nullptr;                \
      }                                                                      \
    }

} // namespace cfdem

#endif // __RUN_TIME_SELECTION_TABLES_H__
