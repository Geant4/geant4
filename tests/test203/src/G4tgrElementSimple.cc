#include "G4tgrElementSimple.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMessenger.hh"
#include "CLHEP/Units/SystemOfUnits.h"

//-------------------------------------------------------------
G4tgrElementSimple::G4tgrElementSimple( const vector<G4String>& wl ) 
{
  //---------- Check for miminum number of words read 
  G4tgrUtils::CheckWLsize( wl, 5, WLSIZE_EQ, "G4tgrElementSimple::G4tgrElementSimple");

  theType = "ElementSimple";
  theName = G4tgrUtils::SubQuotes( wl[1] );
  theSymbol = G4tgrUtils::SubQuotes( wl[2] );
  theZ = G4tgrUtils::GetInt( wl[3] );
  theA = G4tgrUtils::GetFloat( wl[4], g/mole);
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << " G4tgrElementSimple created " << theName << G4endl;
#endif

}

