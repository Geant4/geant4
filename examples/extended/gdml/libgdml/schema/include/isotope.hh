#ifndef ISOTOPE_H
#define ISOTOPE_H 1

#include "MaterialIsotopeType.hh"

class isotope : virtual public SAXObject, public MaterialIsotopeType
{
public:
  isotope() {
  }
  virtual ~isotope() {
  }
  
  virtual SAXObject::Type type() {
    return SAXObject::element;
  }
};

#endif // ISOTOPE_H
