//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4TrackingMessenger.hh 66241 2012-12-13 18:34:42Z gunter $
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

