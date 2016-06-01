// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayXMessenger.hh,v 2.2 1998/07/12 03:44:33 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// GEANT4 Ray tracing for X windows Messenger - John Allison 5th June 1997

#ifdef G4VIS_BUILD_RAYX_DRIVER

#ifndef G4RAYXMESSENGER_HH
#define G4RAYXMESSENGER_HH

#include "G4RayMessenger.hh"

class G4RayX;
class G4UIcommand;

class G4RayXMessenger: public G4RayMessenger {
public:
  G4RayXMessenger (G4RayX* pRayX);
  ~G4RayXMessenger ();
  void SetNewValue (G4UIcommand* command, G4String newValues);
  G4String GetCurrentValue (G4UIcommand* command);
private:
  G4RayX* fpRayX;
};

#endif

#endif
