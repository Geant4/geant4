#ifndef G4ImportanceAlgorithm_hh
#define G4ImportanceAlgorithm_hh 1

#include "G4VImportanceAlgorithm.hh"

class G4ImportanceAlgorithm : public G4VImportanceAlgorithm {

public:
  G4ImportanceAlgorithm();
  ~G4ImportanceAlgorithm();
  G4Nsplit_Weight 
  Calculate(G4double ipre_over_ipost, G4double init_w) const;
private:
  void Error(const G4String &m) const {
    G4cout << "ERROR: G4ImportanceAlgorithm::" << m << G4endl;
    exit(1);
  }
  void Warning(const G4String &m) const {
    G4cout << "WARNING: G4ImportanceAlgorithm::" <<  m << G4endl;
  }
  mutable G4bool fWorned;
};


#endif
