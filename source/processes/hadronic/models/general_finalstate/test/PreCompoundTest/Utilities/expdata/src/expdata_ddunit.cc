#include "expdata_ddunit.h"

ClassImp(expdata_ddunit);



TString expdata_ddunit::GetGraphTitle() const
{
  Char_t theEnergy[5];
  Char_t theAngle[5];
  Char_t target_A[3];
  sprintf(theEnergy,"%.1f",this->GetProjectileEnergy());
  //  sprintf(theAngle,"%.1f",this->GetAngle()*180/3.14159265358979312);
  sprintf(theAngle,"%.1f",this->GetAngle());
  sprintf(target_A,"%d",this->GetTargetA());
  TString Title("Double differential cross section: ^{");
  Title += target_A;
  Title += 
    "}" + this->GetTargetSymbol() + "(" + 
    this->GetProjectileSymbol() + "x)" + 
    this->GetParticleSymbol() + 
    " at " + theEnergy + " MeV " + 
    theAngle + " degrees";
  return Title;
}


void expdata_ddunit::ShowYourSelf(std::ostringstream & os) const 
{
  os << "+-------------+\n"
     << "|    Angle    |  " << angle << '\n'
     << "+-------------+\n\n";
  
  expdata::ShowYourSelf(os);
  return;
}

