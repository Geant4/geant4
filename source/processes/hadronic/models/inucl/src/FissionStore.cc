#include "FissionStore.h"
#include <math.h>

FissionConfiguration FissionStore::generateConfiguration(double amax, 
                    double rand) const {

const double small = -30.;

double totProb = 0.;
vector<double> probs(configurations.size());
//cout << " amax " << amax << " ic " << configurations.size() << endl;
for(int i = 0; i < configurations.size(); i++) {
  double ez = configurations[i].ezet;
  double pr = ez - amax;
  if(pr < small) pr = small;
  pr = exp(pr); 
//  configurations[i].print();
//  cout << " probability " << pr << endl; 
  totProb += pr;
  probs[i] = totProb;  
};

double st = totProb*rand;
int igen = 0;
while (probs[igen] <= st && igen < configurations.size()) igen++;
//cout << " igen " << igen << endl;
return configurations[igen];

}		    
