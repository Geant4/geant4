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

#include "G4UImessenger.hh"
#include "G4OpticalProcessIndex.hh"

#include "globals.hh"

class G4VProcess;
class G4OpticalPhysics;

class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcommand;

// Messenger class that defines commands for the optical physics

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

  void Deprecated(void);

  // data members

  /// associated class
  G4OpticalPhysics*      fOpticalPhysics;

  /// command directory
  G4UIdirectory*         fDir;
  G4UIdirectory*         fDir2;

  /// selected optical process
  G4OpticalProcessIndex  fSelectedProcessIndex;

  /// selectOpProcess command
  G4UIcommand*           fActivateProcessCmd;

  /// setProcessVerbose command
  G4UIcmdWithAnInteger*  fVerboseCmd;

  /// setTrackSecondariesFirst command
  G4UIcommand*           fTrackSecondariesFirstCmd;

  // Cerenkov

  /// setCerenkovMaxPhotons command
  G4UIcmdWithAnInteger*  fCerenkovMaxPhotonsCmd;
  G4UIcmdWithAnInteger*  fCerenkovMaxPhotons1Cmd;

  /// setCerenkovMaxBetaChange command
  G4UIcmdWithADouble*    fCerenkovMaxBetaChangeCmd;
  G4UIcmdWithADouble*    fCerenkovMaxBetaChange1Cmd;

  /// setStackPhotons command
  G4UIcmdWithABool*      fCerenkovStackPhotonsCmd;
  G4UIcmdWithABool*      fCerenkovStackPhotons1Cmd;

  G4UIcmdWithABool*      fCerenkovTrackSecondariesFirstCmd;
  G4UIcmdWithAnInteger*  fCerenkovVerbosityCmd;

  // Scintillation

  /// setScintillationYieldFactor command
  G4UIcmdWithADouble*    fScintYieldFactorCmd;
  G4UIcmdWithADouble*    fScintYieldFactor1Cmd;

  /// setScintillationByParticleType command
  G4UIcmdWithABool*      fScintByParticleTypeCmd;
  G4UIcmdWithABool*      fScintByParticleType1Cmd;

  /// setScintillationTrackInfo command
  G4UIcmdWithABool*      fScintTrackInfoCmd;
  G4UIcmdWithABool*      fScintTrackInfo1Cmd;

  /// setStackPhotons command
  G4UIcmdWithABool*      fScintStackPhotonsCmd;
  G4UIcmdWithABool*      fScintStackPhotons1Cmd;

  G4UIcmdWithADouble*    fScintExcitationRatioCmd;

  G4UIcmdWithABool*      fScintTrackSecondariesFirstCmd;
  G4UIcmdWithAnInteger*  fScintillationVerbosityCmd;

  /// setFiniteRiseTime command
  G4UIcmdWithABool*      fScintFiniteRiseTimeCmd;
  G4UIcmdWithABool*      fScintFiniteRiseTime1Cmd;

  G4UIcmdWithAnInteger*  fScintVerbosityCmd;

  // WLS

  /// setWLSTimeProfile command
  G4UIcmdWithAString*    fWLSTimeProfileCmd;
  G4UIcmdWithAString*    fWLSTimeProfile1Cmd;
  G4UIcmdWithAnInteger*  fWLSVerbosityCmd;

  /// setInvokeSD command
  G4UIcmdWithABool*      fBoundaryInvokeSDCmd;
  G4UIcmdWithABool*      fBoundaryInvokeSD1Cmd;
  G4UIcmdWithAnInteger*  fBoundaryVerbosityCmd;

  G4UIcmdWithAnInteger*  fAbsorptionVerbosityCmd;
  G4UIcmdWithAnInteger*  fRayleighVerbosityCmd;
  G4UIcmdWithAnInteger*  fMieVerbosityCmd;

};

#endif // G4OpticalPhysicsMessenger_h
