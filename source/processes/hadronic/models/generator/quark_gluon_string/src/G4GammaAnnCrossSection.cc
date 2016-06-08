#include "G4GammaAnnCrossSection.hh"
#include "G4Gamma.hh"

G4GammaAnnCrossSection::
G4GammaAnnCrossSection()
{
// pho0 Nucleon
  theGammaNucXSections.push_back(new G4ASCCrossSection(113,  2212, 13.7, 35.9, 0.45, 0.079));
// omega0 Nucleon
  theGammaNucXSections.push_back(new G4ASCCrossSection(223,  2212, 13.7, 35.9, 0.45, 0.079));
// phi0 Nucleon
  theGammaNucXSections.push_back(new G4ASCCrossSection(333,  2212, 12.2, 26.4, 0.50, 0.079));
}

G4bool G4GammaAnnCrossSection::
InCharge(G4int aCode, G4int bCode)
{
  G4bool result = false;
  if(aCode==G4Gamma::Gamma()->GetPDGEncoding())
  {
    result=true;
  }
  else if(bCode==G4Gamma::Gamma()->GetPDGEncoding())
  {
    result = true;
  }
  return result;
}    

G4double G4GammaAnnCrossSection::
GetXsec(G4double s)
{
  G4double result = 0;
  // ratios from Phys.Lett.B40:121-126,1972; 22% assigned to higher resonances

  typedef G4std::vector<G4ASCCrossSection*>::iterator iter;
  iter i;
  for(i=theGammaNucXSections.begin(); i!=theGammaNucXSections.end(); i++)
  {
    result += (*i)->GetXsec(s);
  }
  
  // Account for higher resonances.
  result /= 0.78;

  return result;
}
