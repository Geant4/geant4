#ifndef BOOLEAN_UNION_H
#define BOOLEAN_UNION_H 1

#include "SAXObject.hh"
#include "BooleanSolidType.hh"

class boolean_union : public SAXObject, public BooleanSolidType
{
public:
  boolean_union() {
  }
  virtual ~boolean_union() {
  }
  virtual SAXObject::Type type() {
    return SAXObject::element;
  }
};

#endif // BOOLEAN_UNION_H
