#ifndef FluoTestNormalization_h
#define FluoTestNormalization_h 1
#include "globals.hh"

class FluoTestDataSet;
class FluoTestNormalization
{
public:

  FluoTestNormalization();
  ~FluoTestNormalization();
  const FluoTestDataSet* Normalize(G4double, G4double, G4int,G4String);

private:

  G4double Integrate(G4double, G4double, G4int, FluoTestDataSet*);

};
 
#endif
