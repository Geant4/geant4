#include "InuclSpecialFunctions.h"

double InuclSpecialFunctions::bindingEnergy(double A, double Z) {

// calculates the nuclei binding energy using Kummel or exact or asymptotic
// high temperature 

double DM;
double AN = A - Z;

if(AN < 0.1 || Z < 0.1) {
  DM = 0.;
}
 else {

  if(A <= 256.) {
    if(AN >= 20. && Z >= 20) { 
      if(Z < 1.7*AN && Z > 0.3*AN) { // standard
        DM = bindingEnergyKummel(A,Z);
      }
       else { // bad case
        DM = bindingEnergyAsymptotic(A,Z);
      }; 
    }
     else {
      if(A > 60. || Z > 21) { // bad case
        DM = bindingEnergyAsymptotic(A,Z);
      }
       else { // exact case
        DM = bindingEnergyExact(A,Z);
      }; 
    }; 
  }  
   else {
    DM = bindingEnergyAsymptotic(A,Z);
  }; 
};  

//cout << " A " << A << " Z " << Z << " DM " << DM << endl;
return DM;
}
