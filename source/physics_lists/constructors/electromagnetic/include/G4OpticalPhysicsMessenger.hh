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
//
//---------------------------------------------------------------------------
//
// ClassName:   G4OpticalPhysicsMessenger
//
// Author:      P.Gumplinger 30.09.2009
//
// Modified:    P.Gumplinger 29.09.2011
//              (based on code from I. Hrivnacova)
//
//----------------------------------------------------------------------------
//
// This class defines commands for the optical physics
//

#ifndef G4OpticalPhysicsMessenger_h
#define G4OpticalPhysicsMessenger_h 1

#include "globals.hh"

#include "G4UImessenger.hh"
#include "G4OpticalProcessIndex.hh"

class G4VProcess;
class G4OpticalPhysics;

class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcommand;

// Messenger class that defines commands for the optical physics
//
// Implements commands:
// - /process/optical/processActivation proc_name flag
// - /process/optical/verbose level
// - /process/optical/setTrackSecondariesFirst proc_name flag
// - /process/optical/defaults/cerenkov/setMaxPhotons val
// - /process/optical/defaults/cerenkov/setMaxBetaChange val
// - /process/optical/defaults/scintillation/setYieldFactor val
// - /process/optical/defaults/scintillation/setByParticleType val
// - /process/optical/defaults/scintillation/setFiniteRiseTime val
// - /process/optical/defaults/wls/setTimeProfile val

class G4OpticalPhysicsMessenger: public G4UImessenger
{
  public:

    G4OpticalPhysicsMessenger(G4OpticalPhysics*);
    virtual ~G4OpticalPhysicsMessenger();

    // methods
    virtual void SetNewValue(G4UIcommand*, G4String);
    
private:

  /// Not implemented
  G4OpticalPhysicsMessenger();
  /// Not implemented
  G4OpticalPhysicsMessenger(const G4OpticalPhysicsMessenger& right);
  /// Not implemented
  G4OpticalPhysicsMessenger& operator=(const G4OpticalPhysicsMessenger& right);

  // data members

  /// associated class
  G4OpticalPhysics*     fOpticalPhysics;

  /// command directory
  G4UIdirectory*         fDir;
  G4UIdirectory*         fDir2;

  /// selected optical process
  G4OpticalProcessIndex  fSelectedProcessIndex;

  /// selectOpProcess command
  G4UIcommand*    fActivateProcessCmd;


  /// setProcessVerbose command
  G4UIcmdWithAnInteger*  fSetOpProcessVerboseCmd;

  /// setCerenkovMaxPhotons command
  G4UIcmdWithAnInteger*  fSetCerenkovMaxPhotonsCmd;

  /// setCerenkovMaxBetaChange command
  G4UIcmdWithADouble*    fSetCerenkovMaxBetaChangeCmd;

  /// setScintillationYieldFactor command
  G4UIcmdWithADouble*    fSetScintillationYieldFactorCmd;

  /// setScintillationByParticleType command
  G4UIcmdWithABool*      fSetScintillationByParticleTypeCmd;

//  /// setOpticalSurfaceModel command
//  G4UIcmdWithAString*    fSetOpticalSurfaceModelCmd;

  /// setWLSTimeProfile command
  G4UIcmdWithAString*    fSetWLSTimeProfileCmd;

  /// setTrackSecondariesFirst command
  G4UIcommand*      fSetTrackSecondariesFirstCmd;

  /// setFiniteRiseTime command
  G4UIcmdWithABool*      fSetFiniteRiseTimeCmd;

};

#endif // G4OpticalPhysicsMessenger_h
