// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisManMessLights.cc,v 1.4 1999-12-15 14:54:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// GEANT4 Visualization Manager Messenger - John Allison 22nd July 1996.
// Lights sub-menu.

#include "G4VisManMessenger.hh"

#include "G4VisManager.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4ios.hh"

void G4VisManMessenger::AddCommandLights () {

  G4UIcommand* command;
  G4UIparameter* param;

  ///////////////////////////////////////////////  /vis~/lights/direction ////
  //lights \hline
  //lights /vis~/lights/direction & $\theta$ $\phi$ &
  //lights Set direction from target to light as
  //lights $\theta$, $\phi$ (in degrees). \\%
  command = new G4UIcommand ("/vis~/lights/direction", this);
  command -> SetGuidance
    (
     "Set direction from target to light as "
     "theta, phi (in degrees)."
     );
  param   =  new G4UIparameter ("theta", 'd', true);
  param   -> SetGuidance ("degrees");
  param   -> SetDefaultValue (0.0);
  command -> SetParameter (param);
  param   =  new G4UIparameter ("phi", 'd', true);
  param   -> SetGuidance ("degrees");
  param   -> SetDefaultValue (0.0);
  command -> SetParameter (param);
  fCommandList.append (command);
}

void G4VisManMessenger::DoCommandLights (const G4String& commandPath,
					 G4String& newValues) {

  ////////////////////////////////////////  /vis~/lights/direction  ////
  if (commandPath == "/vis~/lights/direction") {
    G4double theta, phi ;
    const char* aString = newValues;
    G4std::istrstream is((char*) aString) ; is >> theta >> phi;
    theta = theta * deg;
    phi   = phi   * deg;
    G4double x = sin (theta) * cos (phi);
    G4double y = sin (theta) * sin (phi);
    G4double z = cos (theta);
    G4Vector3D lightpointDirection (x, y, z);
    G4ViewParameters& viewParams = fpVMan -> SetCurrentViewParameters();
    viewParams.SetLightpointDirection (lightpointDirection);
    if (fpVMan -> GetVerboseLevel () > 0) {
      G4cout << "Lightpoint direction set to " << lightpointDirection << G4endl;
      if (fpVMan -> GetVerboseLevel () > 1) {
	fpVMan -> PrintCurrentView ();
      }
    }
    if (ViewValid ()) {
      fpVMan -> GetCurrentViewer () -> ClearView ();  // Clears buffer only.
      fpVMan -> Draw ();    
      fpVMan -> Show ();    
    }
  }
}
