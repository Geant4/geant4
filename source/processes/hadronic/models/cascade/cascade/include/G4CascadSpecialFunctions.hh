#include "globals.hh"

#include <math.h>
#include "g4std/algorithm"
#include "g4std/vector"

namespace G4CascadSpecialFunctions {

  G4std::pair<G4int, G4double> getPositionInEnergyScale2(G4double e); 

  G4std::pair<G4int, G4double> getPositionInEnergyScale1(G4double e);
 
  G4double absorptionCrosSection(G4double e, 
				 G4int type);

  G4double crossSection(G4double e, 
			G4int is);

  G4std::pair<G4int, G4double> getPositionInEnergyScaleEMS(G4double e); 

}

