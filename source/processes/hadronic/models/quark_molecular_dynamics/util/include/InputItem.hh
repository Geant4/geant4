#ifndef _InputItem_H
#define _InputItem_H

class istream;
class ostream;
class String;                     
class Boolean;

class InputItem {
  public:
    virtual const String& getKey() const = 0;
    virtual void read(istream& is) = 0;
    virtual void write(ostream& os) = 0;
    virtual Boolean hasBeenSet() const = 0;
    virtual ~InputItem() {}
};

#endif // _InputItem_H





