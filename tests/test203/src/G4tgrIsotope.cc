#include "G4tgrIsotope.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMessenger.hh"
#include "CLHEP/Units/SystemOfUnits.h"

//-------------------------------------------------------------
G4tgrIsotope::G4tgrIsotope( const vector<G4String>& wl ) 
{
  //---------- Check for miminum number of words read 
  G4tgrUtils::CheckWLsize( wl, 5, WLSIZE_EQ, "G4tgrIsotope::G4tgIstotope");

  theName = G4tgrUtils::SubQuotes( wl[1] );
  theZ = G4tgrUtils::GetInt( wl[2] );
  theN = G4tgrUtils::GetInt( wl[3] );
  theA = G4tgrUtils::GetFloat( wl[4], g/mole);
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 1 ) 
    G4cout << " G4tgrIsotope created " << theName << G4endl;
#endif

}

