// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TrackingMessenger.hh,v 1.6 2000-11-11 06:34:10 tsasaki Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4TrackingMessenger.hh
//
// class Description:
//   This is a messenger class to interface to exchange information
//   between tracking/stepping and UI.
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Makoto  Asai   (e-mail: asai@kekvax.kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

#ifndef G4TrackingMessenger_h
#define G4TrackingMessenger_h 1

class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;
class G4TrackingManager;
class G4SteppingManager;
#include "G4UImessenger.hh"
///////////////////////////////////////////////
class G4TrackingMessenger: public G4UImessenger
///////////////////////////////////////////////
{

//--------
public: // without description
//--------

   G4TrackingMessenger(G4TrackingManager* trMan);
   ~G4TrackingMessenger();
   void SetNewValue(G4UIcommand * command,G4String newValues);
   G4String GetCurrentValue(G4UIcommand * command);

//---------
   private:
//---------

   G4TrackingManager* trackingManager;
   G4SteppingManager* steppingManager;

  // commands 
    G4UIdirectory *             TrackingDirectory;
    G4UIcmdWithoutParameter *   AbortCmd;
    G4UIcmdWithoutParameter *   ResumeCmd;
    G4UIcmdWithAnInteger *      StoreTrajectoryCmd;
    G4UIcmdWithAnInteger *      VerboseCmd;

};

#endif

