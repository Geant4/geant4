
#include "G4VStatMFCanonical.hh"

G4VStatMFCanonical & G4VStatMFCanonical::operator=(const G4VStatMFCanonical & right)
{
  G4Exception("G4VStatMFCanonical::operator= meant to not be accessable");
  return *this;
}


G4bool G4VStatMFCanonical::operator==(const G4VStatMFCanonical & right)
{
  return false;
}
 

G4bool G4VStatMFCanonical::operator!=(const G4VStatMFCanonical & right)
{
  return true;
}


G4double G4VStatMFCanonical::Beta(const G4double & T) const 
{
  if (T > G4StatMFParameters::GetCriticalTemp()) return 0.0;
  else {
    G4double CriticalTempSqr = G4StatMFParameters::GetCriticalTemp()*
      G4StatMFParameters::GetCriticalTemp();
    G4double TempSqr = T*T;
    G4double tmp = (CriticalTempSqr-TempSqr)/(CriticalTempSqr+TempSqr);
    return G4StatMFParameters::GetBeta0()*tmp*pow(tmp,1.0/4.0);
  }
}


G4double G4VStatMFCanonical::DBetaDT(const G4double & T) const 
{
  if (T > G4StatMFParameters::GetCriticalTemp()) return 0.0;
  else {
    G4double CriticalTempSqr = G4StatMFParameters::GetCriticalTemp()*
      G4StatMFParameters::GetCriticalTemp();
    G4double TempSqr = T*T;
    G4double tmp = (CriticalTempSqr-TempSqr)/(CriticalTempSqr+TempSqr);
    return -0.5*G4StatMFParameters::GetBeta0()*pow(tmp,1.0/4.0)*
      (CriticalTempSqr*T)/((CriticalTempSqr+TempSqr)*(CriticalTempSqr+TempSqr));
  }
}



void G4VStatMFCanonical::SortFragments(void)
{
  NumOfCharged = 0;

  G4int i;
  for (i = 0; i < Multiplicity; i++) {
    if (FragmentsZ(i) > 0 ) {
      OrderedA.insert(FragmentsA(i));
      OrderedZ.insert(FragmentsZ(i));
      NumOfCharged++;
    }
  }

  NumOfNeutrons = 0;

  for ( i = 0; i < Multiplicity; i++) {
    if (FragmentsZ(i) <= 0 ) {
      OrderedA.insert(FragmentsA(i));
      OrderedZ.insert(0);
      NumOfNeutrons++;
    }
  }

  if (NumOfCharged <= 1) NumOfCharged = Multiplicity;

}
