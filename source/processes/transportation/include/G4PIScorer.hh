#ifndef G4PIScorer_hh
#define G4PIScorer_hh G4PIScorer_hh

#include "G4VPScorer.hh"
#include "G4PMapPtkTallys.hh"

#include "g4std/iostream"
#include "G4PScorer.hh"

class G4Step;
class G4PStep;
class G4VIStore;

class G4PIScorer : public G4VPScorer {
public:
  G4PIScorer(const G4VIStore &IStore);
  ~G4PIScorer();
  void Score(const G4Step &aStep, const G4PStep &aPStep);
  const G4PMapPtkTallys &GetMapPtkTallys() const {
    return fPScorer.GetMapPtkTallys();
  }
  bool CorrectWeight() const {return fCorrectWeight;}
  void ResetCorrectWeight() {fCorrectWeight = true;}
private:
  const G4VIStore &fIStore; 
  G4PScorer fPScorer;
  int fNerr;
  bool fCorrectWeight;
};

G4std::ostream& operator<<(G4std::ostream &out, const G4PIScorer &ps);

#endif



