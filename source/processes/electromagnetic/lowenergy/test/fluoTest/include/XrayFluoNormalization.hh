#ifndef XrayFluoNormalization_h
#define XrayFluoNormalization_h 1
#include "globals.hh"

class XrayFluoDataSet;
class XrayFluoNormalization
{
public:

  XrayFluoNormalization();
  ~XrayFluoNormalization();
  const XrayFluoDataSet* Normalize(G4double, G4double, G4int,G4String);

private:

  G4double Integrate(G4double, G4double, G4int, XrayFluoDataSet*);

};
 
#endif
