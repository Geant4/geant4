#ifndef G4Pscorer_hh
#define G4Pscorer_hh G4Pscorer_hh

#include "G4VPScorer.hh"
#include "G4PMapPtkTallys.hh"

#include <iostream>

class G4Step;
class G4PStep;

class G4PScorer : public G4VPScorer {
public:
  G4PScorer();
  ~G4PScorer();
  void Score(const G4Step &aStep, const G4PStep &aPStep);
  const G4PMapPtkTallys &GetMapPtkTallys() const {
    return fPtkTallys;
  }
private:
  G4PMapPtkTallys fPtkTallys;
};

ostream& operator<<(ostream &out, const G4PScorer &ps);

#endif



