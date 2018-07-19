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
// $Id: DetectorMessenger.hh 103469 2017-04-11 07:29:36Z gcosmo $
//
/// \file medical/GammaTherapy/include/DetectorMessenger.hh
/// \brief Definition of the DetectorMessenger class
//
#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

// -------------------------------------------------------------
//
//
// -------------------------------------------------------------
//      GEANT4 
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- DetectorMessenger -------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of  
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class DetectorMessenger: public G4UImessenger
{
public: // Without description

  DetectorMessenger(DetectorConstruction* );
  virtual ~DetectorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:

  DetectorConstruction* fDetector;
    
  G4UIdirectory*             fDetDir;
  G4UIdirectory*             fDetDir1;
  G4UIdirectory*             fDetDir2;

  G4UIcmdWithAString*        fAbsMaterCmd;
  G4UIcmdWithADoubleAndUnit* fAbsThickCmd;
  G4UIcmdWithADoubleAndUnit* fAbsGapCmd;
  G4UIcmdWithADoubleAndUnit* fAbsSizYZCmd;
  G4UIcmdWithAString*        fWorldMaterCmd;
  G4UIcmdWithADoubleAndUnit* fWorldXCmd;
  G4UIcmdWithADoubleAndUnit* fXMagFieldCmd;
  G4UIcmdWithADoubleAndUnit* fYMagFieldCmd;
  G4UIcmdWithADoubleAndUnit* fZMagFieldCmd;
  G4UIcmdWithAnInteger*      fNumOfAbsCmd;
  G4UIcmdWithAnInteger*      fNumOfEvt;
  G4UIcmdWithAnInteger*      fVerbCmd;
  G4UIcmdWithAnInteger*      fIntCmd;
  G4UIcmdWithADoubleAndUnit* fDeltaECmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif





