#ifndef B01scorer_hh
#define B01scorer_hh B01scorer_hh

#include "G4VPScorer.hh"
#include "G4PMapPtkTallys.hh"

#include <iostream>

class G4Step;
class G4PStep;

class B01Scorer : public G4VPScorer {
public:
  B01Scorer();
  ~B01Scorer();
  void Score(const G4Step &aStep, const G4PStep &aPStep);
  const G4PMapPtkTallys &GetMapPtkTallys() const {
    return fPtkTallys;
  }
private:
  G4PMapPtkTallys fPtkTallys;
};

ostream& operator<<(ostream &out, const B01Scorer &ps);

#endif



