#ifndef DORDREF_H
#define DORDREF_H 1

#include "QuantityType.hh"
#include "ReferenceType.hh"
#include "TagorTagref.hh"

class D : public SAXObject, public QuantityType
{
public:
  D() {
    set_type( "Density" );
    set_unit( "g/mole" );
  }
  virtual ~D() {
  }
  virtual SAXObject::Type type() {
    return SAXObject::element;
  }
};

class Dref : public SAXObject, public ReferenceType
{
public:
  Dref() {
  }
  virtual ~Dref() {
  }
  virtual SAXObject::Type type() {
    return SAXObject::element;
  }
};

typedef TagorTagref DorDref;

#endif // DORDREF_H
