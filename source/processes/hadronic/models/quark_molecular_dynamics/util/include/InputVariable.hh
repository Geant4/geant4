#ifndef _InputVariable_H
#define _InputVariable_H

#ifdef IS_GCC
#pragma interface
#endif

#include "g4std/iostream"

#include "InputItem.hh"
#include "Fallible.hh"
#include "String.hh"

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
    virtual void read(G4std::istream& is) {
      T value;
      is >> value;
      if( is.good() ) {
	validate(value);
      } 
    }
    virtual void write(G4std::ostream& os) { 
      if( isValid() ) {
	os << T(*this);
      }
    }
    virtual Boolean hasBeenSet() const { return isValid(); }            
    
  private: 
    const String theKey;
};

#endif // _InputVariable_H

