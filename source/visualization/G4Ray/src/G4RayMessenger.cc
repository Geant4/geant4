// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayMessenger.cc,v 2.4 1998/07/13 17:11:10 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// GEANT4 Ray tracing for X windows Messenger - John Allison 5th June 1997

#include "G4RayMessenger.hh"

#include "G4RayScene.hh"
#include "G4RayView.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"

G4RayMessenger::G4RayMessenger ():
  fpScene (0),
  fpView (0)
{
  G4UIcommand* command;
  G4UIparameter* param;

  /////////////////////////////////////////////////  /vis~/ray/...  ////
  //vis \hline
  //vis /vis~/ray/ &&
  //vis ...menu of ray tracing commands. \\%
  command = new G4UIcommand ("/vis~/ray/", this);
  command -> SetGuidance ("GEANT4 Ray Tracing");
  fpVisRayCommandDirectory = command;

  /////////////////////////////////////////////////  /vis~/ray/resolution  ////
  //vis \hline
  //vis /vis~/ray/resolution/ & resolution &
  //vis Set resolution in pixels.  A resolution less than 1 causes
  //vis over-sampling which improves the appearance and avoids random
  //vis spots and staircase effects (aliasing). \\%
  command = new G4UIcommand ("/vis~/ray/resolution", this);
  command -> SetGuidance
    ("Set resolution in pixels.  A resolution less than 1 causes"
     "\nover-sampling which improves the appearance and avoids random"
     "\nspots and staircase effects (aliasing)."
     );
  param   =  new G4UIparameter ("resolution", 'd', true);
  param   -> SetDefaultValue  (-1);
  command -> SetParameter     (param);
  fpResolutionCommand = command;
}

G4RayMessenger::~G4RayMessenger () {
  delete fpResolutionCommand;
  delete fpVisRayCommandDirectory;
}

void G4RayMessenger::SetNewValue (G4UIcommand* command, G4String newValues)
{
  /////////////////////////////////////////////////  /vis~/ray/resolution  ////
  if (command == fpResolutionCommand) {
    G4double resolution;
    const char* t = newValues;
    istrstream is ((char*)t); is >> resolution;
    if (resolution < 0) {
      G4cout << "Resolution is currently "
	   << fpView -> GetResolution ()
	   << "\nSupply required resolution as argument."
	   << endl;
    }
    else {
      fpView -> SetResolution (resolution);
      G4cout << "Resolution is now "
	   << fpView -> GetResolution ()
	   << endl;
    }
  }
}

G4String G4RayMessenger::GetCurrentValue (G4UIcommand* command) {
  return "";
}
