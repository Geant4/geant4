//
// File name:     RadmonDetectorMessenger.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMessenger.cc,v 1.1 2005-09-09 08:26:24 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorMessenger.hh"

#include "globals.hh"



                                                RadmonDetectorMessenger :: RadmonDetectorMessenger(RadmonVDetectorLayout * layout)
:
 detectorLayout(layout)
{
 // TO BE DONE
 G4cout << "RadmonDetectorMessenger::RadmonDetectorMessenger(): NOT IMPLEMENTED YED" << G4endl;
}



                                                RadmonDetectorMessenger :: ~RadmonDetectorMessenger()
{
 // TO BE DONE
 G4cout << "RadmonDetectorMessenger::~RadmonDetectorMessenger(): NOT IMPLEMENTED YED" << G4endl;
}





G4String                                        RadmonDetectorMessenger :: GetCurrentValue(G4UIcommand * /* command */)
{
 // TO BE DONE
 G4cout << "RadmonDetectorMessenger::GetCurrentValue(): NOT IMPLEMENTED YED" << G4endl;
 
 return G4String();
}



void                                            RadmonDetectorMessenger :: SetNewValue(G4UIcommand * /* command */, G4String /* newValue */)
{
 // TO BE DONE
 G4cout << "RadmonDetectorMessenger::SetNewValue(): NOT IMPLEMENTED YED" << G4endl;
}
