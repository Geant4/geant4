#ifndef _InputItem_H
#define _InputItem_H

#include "globals.hh"
#include "g4std/iostream"

class String;                     
class Boolean;

class InputItem {
  public:
    virtual const String& getKey() const = 0;
    virtual void read(G4std::istream& is) = 0;
    virtual void write(G4std::ostream& os) = 0;
    virtual Boolean hasBeenSet() const = 0;
    virtual ~InputItem() {}
};

#endif // _InputItem_H





