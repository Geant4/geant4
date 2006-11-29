#include "G4tgbIsotope.hh"
#include "globals.hh"
#include "G4tgrMessenger.hh"

//----------------------------------------------------------------------
G4tgbIsotope::G4tgbIsotope( G4tgrIsotope* hg )
{
  theTgrIsot = hg;
  theG4Isot = 0;
}


//----------------------------------------------------------------------
G4Isotope* G4tgbIsotope::BuildG4Isotope()
{
  G4Isotope* isot = 0;

  //-------- if G4Isotope not found, construct it 
  if( theG4Isot == 0 ) { 
  //----- construct new G4Isotope 
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " constructing  new G4Isotope " <<  " " << theTgrIsot->GetName() << " " << theTgrIsot->GetZ() << " " << theTgrIsot->GetN() << " " << theTgrIsot->GetA() << endl;
#endif
    isot = new G4Isotope(theTgrIsot->GetName(), theTgrIsot->GetZ(), theTgrIsot->GetN(), theTgrIsot->GetA() );
    theG4Isot = isot;
  } else {
    isot = theG4Isot; 
  }

  return isot;
}

