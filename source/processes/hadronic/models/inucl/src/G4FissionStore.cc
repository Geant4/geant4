#include "G4FissionStore.h"
#include <math.h>

G4FissionConfiguration G4FissionStore::generateConfiguration(G4double amax, 
                    G4double rand) const {

const G4double small = -30.0;

G4double totProb = 0.0;
vector<G4double> probs(configurations.size());
// G4cout << " amax " << amax << " ic " << configurations.size() << G4endl;
for(G4int i = 0; i < configurations.size(); i++) {
  G4double ez = configurations[i].ezet;
  G4double pr = ez - amax;
  if(pr < small) pr = small;
  pr = exp(pr); 
//  configurations[i].print();
//  G4cout << " probability " << pr << G4endl; 
  totProb += pr;
  probs[i] = totProb;  
};

G4double st = totProb * rand;
G4int igen = 0;
while (probs[igen] <= st && igen < configurations.size()) igen++;
// G4cout << " igen " << igen << G4endl;
return configurations[igen];

}		    
