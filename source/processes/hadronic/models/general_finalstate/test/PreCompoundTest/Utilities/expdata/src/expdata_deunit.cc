#include "expdata_deunit.h"
#include <stdlib.h>

ClassImp(expdata_deunit);


TString expdata_deunit::GetGraphTitle() const 
{
  Char_t theEnergy[10];
  sprintf(theEnergy,"%.1f",this->GetProjectileEnergy());
  Char_t target_A[3];
  sprintf(target_A,"%d",Int_t(this->GetTargetA()));
  TString Title("Differential Cross Section: ^{");
  Title += target_A;
  Title +=
    "}" + this->GetTargetSymbol() + 
    "(" + this->GetProjectileSymbol() + "," + "x)" + this->GetParticleSymbol() +
    " at " + theEnergy + " MeV ";
  return Title;
}


void expdata_deunit::ShowYourSelf(std::ostringstream & os) const 
{
  os << "****************************************************************\n"
     << "*                   ENERGY DISTRIBUTION                        *\n" 
     << "****************************************************************\n\n";

  expdata::ShowYourSelf(os);
  os << "\n\tThere are " << this->GetN() << " points\n";
  return;
}
