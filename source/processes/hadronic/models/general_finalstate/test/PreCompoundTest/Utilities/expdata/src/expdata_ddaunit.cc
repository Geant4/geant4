#include "expdata_ddaunit.h"

ClassImp(expdata_ddaunit);

TString expdata_ddaunit::GetGraphTitle() const
{
  Char_t theEnergyMin[5];
  Char_t theEnergyMax[5];
  Char_t theEnergy[5];
  Char_t theTargetA[5];
  sprintf(theEnergy,"%.1f",this->GetProjectileEnergy());
  sprintf(theEnergyMin,"%.1f",this->E_min);
  sprintf(theEnergyMax,"%.1f",this->E_max);
  sprintf(theTargetA,"%d",this->GetTargetA());
  TString Title("Differential cross section: ^{");
  Title += theTargetA;
  Title += 
    "}" + this->GetTargetSymbol() + "(" +
    this->GetProjectileSymbol() + "x)" +
    this->GetParticleSymbol() +
    " at " + theEnergy + " MeV  #DeltaE = [" +
    theEnergyMin + ',' + theEnergyMax + "] MeV";
  return Title;
}

void expdata_ddaunit::ShowYourSelf(std::ostringstream & os) const
{
  os << "+-------------+\n"
     << "| Energy range|  [" << E_min << ',' << E_max << "] MeV\n"
     << "+-------------+\n\n";
  
  expdata::ShowYourSelf(os);
  return;
}
