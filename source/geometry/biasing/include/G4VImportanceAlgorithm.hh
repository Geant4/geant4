#ifndef G4VImportanceAlgorithm_hh
#define G4VImportanceAlgorithm_hh 1

#include "G4Nsplit_Weight.hh"

class G4VImportanceAlgorithm {
public:
  virtual ~G4VImportanceAlgorithm(){}
  virtual G4Nsplit_Weight 
  Calculate(G4double ipre_over_ipost, G4double init_w) const = 0;
};



#endif

