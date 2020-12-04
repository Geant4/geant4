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

#ifndef G4OpticalParametersMessenger_h
#define G4OpticalParametersMessenger_h 1

#include "G4UImessenger.hh"
//#include "G4OpticalProcessIndex.hh"

#include "globals.hh"

class G4VProcess;
class G4OpticalParameters;

class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcommand;

// Messenger class that defines commands for the optical physics

class G4OpticalParametersMessenger: public G4UImessenger
{
  public:

    G4OpticalParametersMessenger(G4OpticalParameters*);
    virtual ~G4OpticalParametersMessenger();

    // methods
    virtual void SetNewValue(G4UIcommand*, G4String);
    
private:

  G4OpticalParametersMessenger() = delete;
  G4OpticalParametersMessenger(const G4OpticalParametersMessenger& right) = delete;
  G4OpticalParametersMessenger& operator=(const G4OpticalParametersMessenger& right) = delete;

  void Deprecated(void);

  // data members

  /// associated class
  G4OpticalParameters*      params;

  /// command directory
  G4UIdirectory*         fDir;
  G4UIdirectory*         fDir2;

  /// selected optical process
  //G4OpticalProcessIndex  fSelectedProcessIndex;

  /// selectOpProcess command
  G4UIcommand*           fActivateProcessCmd;

  /// setProcessVerbose command
  G4UIcmdWithAnInteger*  fVerboseCmd;

  /// setTrackSecondariesFirst command
  G4UIcommand*           fTrackSecondariesFirstCmd;

  // Cerenkov

  // setCerenkovMaxPhotons command
  G4UIcmdWithAnInteger*  fCerenkovMaxPhotonsCmd;
  G4UIcmdWithAnInteger*  fCerenkovMaxPhotons1Cmd;

  /// setCerenkovMaxBetaChange command
  G4UIcmdWithADouble*    fCerenkovMaxBetaChangeCmd;
  G4UIcmdWithADouble*    fCerenkovMaxBetaChange1Cmd;

  /// setStackPhotons command
  G4UIcmdWithABool*      fCerenkovStackPhotonsCmd;
  G4UIcmdWithABool*      fCerenkovStackPhotons1Cmd;

  G4UIcmdWithABool*      fCerenkovTrackSecondariesFirstCmd;
  G4UIcmdWithAnInteger*  fCerenkovVerboseLevelCmd;

  // Scintillation

  // setScintillationYieldFactor command
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

  /// setFiniteRiseTime command
  G4UIcmdWithABool*      fScintFiniteRiseTimeCmd;
  G4UIcmdWithABool*      fScintFiniteRiseTime1Cmd;

  // setEnhancedTimeConstants commnad
  G4UIcmdWithABool*      fScintEnhancedTimeConstantsCmd;

  G4UIcmdWithAnInteger*  fScintVerboseLevelCmd;

  // WLS

  /// setWLSTimeProfile command
  G4UIcmdWithAString*    fWLSTimeProfileCmd;
  G4UIcmdWithAString*    fWLSTimeProfile1Cmd;
  G4UIcmdWithAnInteger*  fWLSVerboseLevelCmd;

  // WLS2

  /// setWLS2TimeProfile command
  G4UIcmdWithAString*    fWLS2TimeProfileCmd;
  G4UIcmdWithAnInteger*  fWLS2VerboseLevelCmd;
  
  /// setInvokeSD command
  G4UIcmdWithABool*      fBoundaryInvokeSDCmd;
  G4UIcmdWithABool*      fBoundaryInvokeSD1Cmd;
  G4UIcmdWithAnInteger*  fBoundaryVerboseLevelCmd;

  G4UIcmdWithAnInteger*  fAbsorptionVerboseLevelCmd;
  G4UIcmdWithAnInteger*  fRayleighVerboseLevelCmd;
  G4UIcmdWithAnInteger*  fMieVerboseLevelCmd;

  G4UIcommand*           fDumpCmd;
};

#endif // G4OpticalParametersMessenger_h
