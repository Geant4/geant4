#include "G4tgrElementFromIsotopes.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMessenger.hh"
#include "CLHEP/Units/SystemOfUnits.h"

//-------------------------------------------------------------
G4tgrElementFromIsotopes::G4tgrElementFromIsotopes( const vector<G4String>& wl ) 
{
  //---------- Check for miminum number of words read 
  G4tgrUtils::CheckWLsize( wl, 6, WLSIZE_LE, "G4tgrElementFromIsotopes::G4tgrElementFromIsotopes");
    //:ELEM_FROM_ISOT NAME SYMBOL N_ISOT (ISOT_NAME ISOT_ABUNDANCE)

  theType = "ElementFromIsotopes";
  theName = G4tgrUtils::SubQuotes( wl[1] );
  theSymbol = G4tgrUtils::SubQuotes( wl[2] );
  theNoIsotopes = G4tgrUtils::GetInt( wl[3] );

  for( size_t ii = 0; ii < theNoIsotopes; ii++ ){
    theComponents.push_back( G4tgrUtils::SubQuotes( wl[4+ii*2] ) ); 
    theAbundances.push_back( G4tgrUtils::GetFloat( wl[4+ii*2+1] ) );
  }

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
    G4cout << " G4tgrElementFromIsotopes created " << theName << G4endl;
#endif

}

