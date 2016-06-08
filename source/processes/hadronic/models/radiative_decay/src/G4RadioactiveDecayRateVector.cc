

#include "G4RadioactiveDecayRateVector.hh"


G4RadioactiveDecayRateVector::G4RadioactiveDecayRateVector()
{
  ;
  //do nothing at the momment
}



G4RadioactiveDecayRateVector::G4RadioactiveDecayRateVector(const G4RadioactiveDecayRateVector &right)
{
  ionName = right.ionName;
  itsRates = right.itsRates;
}

G4RadioactiveDecayRateVector & G4RadioactiveDecayRateVector::operator=(const G4RadioactiveDecayRateVector &right)
{
  if (this != &right) { 
    ionName = right.ionName;
    itsRates = right.itsRates;
  }
  return *this;
}


G4RadioactiveDecayRateVector::~G4RadioactiveDecayRateVector()
{ ;} 










