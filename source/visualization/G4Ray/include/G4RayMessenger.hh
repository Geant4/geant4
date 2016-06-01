// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayMessenger.hh,v 2.2 1998/07/12 03:44:31 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// GEANT4 Ray tracing Messenger - John Allison 5th June 1997

#ifndef G4RAYMESSENGER_HH
#define G4RAYMESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"

class G4RayX;
class G4RayScene;
class G4RayView;

class G4UIcommand;

class G4RayMessenger: public G4UImessenger {
public:
  G4RayMessenger ();
  ~G4RayMessenger ();
  void SetNewValue (G4UIcommand* command, G4String newValues);
  G4String GetCurrentValue (G4UIcommand* command);
  void RegisterScene (G4RayScene* pScene);
  void RegisterView (G4RayView* pView);
private:
  G4RayScene* fpScene;
  G4RayView* fpView;
  G4UIcommand* fpVisRayCommandDirectory;
  G4UIcommand* fpResolutionCommand;
};

#include "G4RayMessenger.icc"

#endif
