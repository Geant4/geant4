#ifndef B02scorer_hh
#define B02scorer_hh B02scorer_hh

#include "G4VPScorer.hh"
#include "G4PMapPtkTallys.hh"

#include <iostream>

class G4Step;
class G4PStep;

class B02Scorer : public G4VPScorer {
public:
  B02Scorer();
  ~B02Scorer();
  void Score(const G4Step &aStep, const G4PStep &aPStep);
  const G4PMapPtkTallys &GetMapPtkTallys() const {
    return fPtkTallys;
  }
private:
  G4PMapPtkTallys fPtkTallys;
};

ostream& operator<<(ostream &out, const B02Scorer &ps);

#endif



