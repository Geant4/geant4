#ifndef ATOM_HH
#define ATOM_HH 1

#include "SAXObject.hh"
#include "AtomType.hh"

class atom : public SAXObject, public AtomType
{
public:
  atom() {
  }
  virtual ~atom() {
  }
  virtual SAXObject::Type type() {
    return SAXObject::element;
  }
};

#endif // ATOM_HH
