#include "G4tgrMaterialMixture.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMessenger.hh"

#include "CLHEP/Units/SystemOfUnits.h"

//-----------------------------------------------------------
G4tgrMaterialMixture::G4tgrMaterialMixture(const G4String& matType, const vector<G4String>& wl)
{
  //---------- Check for miminum number of words read 
  G4tgrUtils::CheckWLsize( wl, 6, WLSIZE_GE, "G4tgrMaterialMixture::G4tgrMaterialMixture");

  theMateType = matType;
  
  //---------- Fill private data 
  theName = G4tgrUtils::SubQuotes( wl[1] );
  theDensity = fabs(G4tgrUtils::GetFloat( wl[2], g/cm3 ) );
  theNoComponents = G4tgrUtils::GetInt( wl[3] );

  G4tgrUtils::CheckWLsize( wl, 4+theNoComponents*2, WLSIZE_GE, "G4tgrMaterialMixture::G4tgrMaterialMixture");
  for(int ii=0; ii<theNoComponents; ii++){
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
         G4cout << ii << " G4tgrMaterialMixture adding component " << wl[ii*2+4] << G4endl;
#endif
    theComponents.push_back(  G4tgrUtils::SubQuotes( wl[ii*2+4] ) );
    theFractions.push_back( G4tgrUtils::GetFloat(wl[ii*2+1+4]) );
  }

  //--------- transform fractions to fractions by weight;
  //-  TransformToFractionsByWeight();
}


//-----------------------------------------------------------
ostream& operator<<(ostream& os, const G4tgrMaterialMixture& mate) 
{
  os << "MATERIAL: " << mate.theName << G4endl
     << "density= " << mate.theDensity << " g/cm3. Number of Components: "
     << mate.theNoComponents << G4endl;
  for (int ii=0; ii<mate.theNoComponents; ii++)
    os << '\t' << mate.theComponents[ii] << '\t' << mate.theFractions[ii] << G4endl;
  return os;
}
