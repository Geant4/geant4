#ifndef BOOLEAN_INTERSECTION_H
#define BOOLEAN_INTERSECTION_H 1

#include "SAXObject.hh"
#include "BooleanSolidType.hh"

class boolean_intersection : public SAXObject, public BooleanSolidType
{
public:
  boolean_intersection() {
  }
  virtual ~boolean_intersection() {
  }
  virtual SAXObject::Type type() {
    return SAXObject::element;
  }
};

#endif // BOOLEAN_INTERSECTION_H
