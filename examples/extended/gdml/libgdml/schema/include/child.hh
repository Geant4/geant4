#ifndef CHILD_H
#define CHILD_H 1

#include "SAXObject.hh"

#include "SinglePlacementType.hh"

class child : public SAXObject, public SinglePlacementType
{
public:
  child() {
  }
  virtual ~child() {
  }
  virtual SAXObject::Type type() {
    return SAXObject::element;
  }
};



#endif // CHILD_H
