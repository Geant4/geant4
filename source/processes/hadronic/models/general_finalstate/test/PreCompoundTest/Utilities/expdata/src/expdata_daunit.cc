#include "expdata_daunit.h"
#include <stdlib.h>

ClassImp(expdata_daunit);


TString expdata_daunit::GetGraphTitle() const 
{
  Char_t theEnergy[10];
  sprintf(theEnergy,"%.1f",this->GetProjectileEnergy());
  Char_t target_A[3];
  sprintf(target_A,"%d",this->GetTargetA());
  TString Title("Differential Cross Section: ^{");
  Title += target_A;
  Title += 
    "}" + this->GetTargetSymbol() + 
    "(" + this->GetProjectileSymbol() + "," + "x)" + this->GetParticleSymbol() +
    " at " + theEnergy + " MeV ";
  return Title;
}


void expdata_daunit::ShowYourSelf(std::ostringstream & os) const 
{
  os << "****************************************************************\n"
     << "*                  ANGULAR DISTRIBUTION                        *\n" 
     << "****************************************************************\n\n";

  expdata::ShowYourSelf(os);
  os << "\n\tThere are " << this->GetN() << " points\n";
  return;
}
