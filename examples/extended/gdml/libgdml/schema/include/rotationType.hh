#ifndef ROTATIONTYPE_H
#define ROTATIONTYPE_H 1


#include "IdentifiableQuantityVectorType.hh"


class rotationType : public IdentifiableQuantityVectorType {
public:
  rotationType() {
    set_unit( "radian" );
    set_type( "cartesian" );
  }
  ~rotationType() {
  }
};



#endif // ROTATIONTYPE_H
