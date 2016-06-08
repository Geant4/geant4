#ifndef rw2stl_algo_h
#define rw2stl_algo_h

#include "rw/defs.h"

template <class T>
class RW2STL_LessPtr
{
public:

  RWBoolean operator()(const T* a, const T* b) const
    { 
      if(a==0) return false;
      if(b==0) return true;
      return *a<*b;
    }

};

#endif
