#ifndef POSITIONTYPE_H
#define POSITIONTYPE_H 1


#include "IdentifiableQuantityVectorType.hh"


class positionType : public IdentifiableQuantityVectorType {
public:
  positionType() {
    set_unit( "mm" );
    set_type( "cartesian" );
  }
  ~positionType() {
  }
};



#endif // POSITIONTYPE_H
