#ifndef ELEMENT_H
#define ELEMENT_H 1

#include "MaterialElementType.hh"

class element : virtual public SAXObject, public MaterialElementType
{
public:
  element() {
  }
  virtual ~element() {
  }
  
  virtual SAXObject::Type type() {
    return SAXObject::element;
  }
};

#endif // ELEMENT_H
