#include <math.h>
#include "pair.h"
#include "vector"

namespace G4InuclSpecialFunctions {

G4double bindingEnergyExact(G4double A, G4double Z);

G4double bindingEnergyKummel(G4double A, G4double Z);

G4double bindingEnergy(G4double A, G4double Z);

G4double bindingEnergyAsymptotic(G4double A, G4double Z);

G4double FermiEnergy(G4double A, G4double Z, G4int ntype);
  
pair<vector<G4double>, vector<G4double> > paraMaker(G4double Z);

pair<G4double, G4double> paraMakerTruncated(G4double Z); 

G4double getAL(G4double A);
 
G4double csNN(G4double e);

G4double csPN(G4double e);

G4double inuclRndm();

G4double randomGauss(G4double sigma);

G4double randomPHI();

pair<G4double, G4double> randomCOS_SIN();

G4double nucleiLevelDensity(G4double a);

vector<G4double> generateWithFixedTheta(G4double ct, G4double p);

}
