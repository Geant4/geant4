// Rich advanced example for Geant4
// AerogelTypeSpec.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef AerogelTypeSpec_h
#define AerogelTypeSpec_h 1
//  different types of Aerogels.
static const G4int MaxNumberOfAerogelTiles=5;
static const G4int MaxNumberOfAerogelTypes=5;
enum AerogelType{AerogelTypeA,AerogelTypeB,AerogelTypeC,AerogelTypeD,
                AerogelTypeE};
static const G4String AerogelTypeString[]={"AerogelTypeA","AerogelTypeB",
					   "AerogelTypeC","AerogelTypeD",
                                           "AerogelTypeE"};
// Type A has dimensions 7*8*4 cm
// In the G4Example only 1 type is used. In the LHCb implementation
// there are 5 types of aerogel.
// In the G4Example only 1 tile is implemented. In the LHCb
// implementation upto 5 tiles are allowed.
#endif 
