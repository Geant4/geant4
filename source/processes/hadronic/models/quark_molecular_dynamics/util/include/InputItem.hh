#ifndef _InputItem_H
#define _InputItem_H

#include "globals.hh"
#include <iostream>

class String;                     
class Boolean;

class InputItem {
  public:
    virtual const String& getKey() const = 0;
    virtual void read(std::istream& is) = 0;
    virtual void write(std::ostream& os) = 0;
    virtual Boolean hasBeenSet() const = 0;
    virtual ~InputItem() {}
};

#endif // _InputItem_H





