#include <math.h>

#include "pair.h"
#include "vector"

namespace G4CascadSpecialFunctions {

  pair<G4int, G4double> getPositionInEnergyScale2(G4double e); 

  pair<G4int, G4double> getPositionInEnergyScale1(G4double e);
 
  G4double absorptionCrosSection(G4double e, 
				 G4int type);

  G4double crossSection(G4double e, 
			G4int is);

  pair<G4int, G4double> getPositionInEnergyScaleEMS(G4double e); 
 
}

