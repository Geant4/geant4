#ifndef ATOMTYPE_H
#define ATOMTYPE_H 1


#include "QuantityType.hh"


class AtomType : public QuantityType {
public:
  AtomType() {
    set_unit("g/mole");
    set_type("A");
  }
  ~AtomType() {
  }
};



#endif // ATOMTYPE_H
