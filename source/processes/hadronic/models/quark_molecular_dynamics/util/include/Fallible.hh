#ifndef _Fallible_H
#define _Fallible_H

#ifdef IS_GCC
#pragma interface
#endif

#include "Definitions.hh"
#include "FallibleBase.hh"

template<class T>
class Fallible
  : private FallibleBase
{
  public:
    Fallible(const T& D) : FallibleBase(Boolean::False),
                   DefaultValue(D),DefaultSet(Boolean::True) {}
    Fallible() : FallibleBase(Boolean::False),DefaultSet(Boolean::False) {}

    FallibleBase::isInvalid;
    FallibleBase::isValid;
    FallibleBase::invalidate;
    T& validate(const T& t) {instance = t; FallibleBase::validate(); 
                             return instance; }
    
    operator T() const {
      if( isInvalid() && !DefaultSet ) {
	FallibleBase::throwError();
      }
      return isValid() ? instance : DefaultValue;
    }

    
    T elseDefaultTo(const T& defaultValue) const { 
      return isValid() ? instance : defaultValue;
    }

    
  protected:
    T getValue() { return isValid() ? instance : DefaultValue; }
  private:
    T instance;  
    T DefaultValue;
    Boolean DefaultSet;
};

#endif // _Fallible_H





