#include "G4tgrMaterialSimple.hh"
#include "G4tgrUtils.hh"

#include "CLHEP/Units/SystemOfUnits.h"


//-------------------------------------------------------------
G4tgrMaterialSimple::G4tgrMaterialSimple(const G4String& matType, const vector<G4String>& wl)
{
  //---------- Check for miminum number of words read 
  G4tgrUtils::CheckWLsize( wl, 5, WLSIZE_EQ, "G4tgrMaterialSimple::G4tgrMaterialSimple");

  theMateType = matType;

  //---------- Fill data from 'wl'
  //-  G4tgrUtils::DumpVS(wl, "G4tgrMaterialSimple::FillData");

  //---------- Fill private data 
  theName = G4tgrUtils::SubQuotes( wl[1] );
  theZ = G4tgrUtils::GetFloat( wl[2] );
  theA = G4tgrUtils::GetFloat( wl[3], g/mole);
  theDensity = G4tgrUtils::GetFloat( wl[4], g/cm3);
  theNoComponents = 0;

}


//-------------------------------------------------------------
ostream& operator<<(ostream& os, const G4tgrMaterialSimple& mate) 
{
  os << "MATERIAL SIMPLE: " << mate.theName 
     << " Z " << mate.theZ << " A " << mate.theA 
     << "density= " << mate.theDensity*g/cm3 << " g/cm3. Number of Components: "
     << mate.theNoComponents << G4endl;
  return os;
}


