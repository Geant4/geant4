#ifndef ELEMENT_H
#define ELEMENT_H 1

#include "MaterialMixtureType.hh"

class material : virtual public SAXObject, public MaterialMixtureType
{
public:
  material() {
  }
  virtual ~material() {
  }
  
  virtual SAXObject::Type type() {
    return SAXObject::element;
  }
};

#endif // ELEMENT_H
