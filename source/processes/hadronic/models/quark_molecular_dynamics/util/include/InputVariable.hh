#ifndef _InputVariable_H
#define _InputVariable_H

#ifdef IS_GCC
#pragma interface
#endif

#include <iostream.h>

#include "InputItem.hh"
#include "Fallible.hh"
#include "String.hh"

class istream;

template<class T>
class InputVariable
  : public InputItem,
    public Fallible<T>
{ 
  public:
    InputVariable(const String key) : theKey(key) {}
    InputVariable(const String key,const T& Default) 
      : theKey(key),Fallible<T>(Default) {}
    
    virtual const String& getKey() const { return theKey; }
    virtual void read(istream& is) {
      T value;
      is >> value;
      if( is.good() ) {
	validate(value);
      } 
    }
    virtual void write(ostream& os) { 
      if( isValid() ) {
	os << T(*this);
      }
    }
    virtual Boolean hasBeenSet() const { return isValid(); }            
    
  private: 
    const String theKey;
};

#endif // _InputVariable_H

