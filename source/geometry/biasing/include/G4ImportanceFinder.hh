#ifndef G4ImportanceFinder_hh
#define G4ImportanceFinder_hh G4ImportanceFinder_hh 

#include "globals.hh"

class G4VParallelStepper;
class G4VIStore;
class G4PTouchableKey;

class G4ImportanceFinder{
public:
  G4ImportanceFinder(const G4VIStore &aIStore):
    fIStore(aIStore)
  {}
  ~G4ImportanceFinder(){}
  
  G4double GetIPre_over_IPost(const G4PTouchableKey &prekey,
			      const G4PTouchableKey &postkey)  const;
private:
  void Error(const G4String &m) const {
    G4cout << "G4ImportanceFinder::" << m << G4endl;
    exit(1);
  }
  const G4VIStore &fIStore;
};

#endif




