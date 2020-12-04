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
/// \file SAXSDetectorConstructionMessenger.hh
/// \brief Implementation of the SAXSDetectorConstructionMessenger class
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef SAXSDetectorConstructionMessenger_h
#define SAXSDetectorConstructionMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class SAXSDetectorConstruction;

class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// DetectorConstruction messenger.

class SAXSDetectorConstructionMessenger: public G4UImessenger
{
public:
    SAXSDetectorConstructionMessenger(SAXSDetectorConstruction* detconstr);
    ~SAXSDetectorConstructionMessenger();

    virtual void SetNewValue(G4UIcommand* command,G4String newValues);

private:
    SAXSDetectorConstruction* fDetector;
    
    G4UIdirectory* fCmdDir;   
    
        G4UIcmdWithAString*        fSetCustomMatFFfilename;        
        G4UIcmdWithADouble* fSetCustomMatDensityCmd;
        G4UIcmdWithADouble* fSetCustomMatHmassfractCmd;
        G4UIcmdWithADouble* fSetCustomMatCmassfractCmd;
        G4UIcmdWithADouble* fSetCustomMatNmassfractCmd;
        G4UIcmdWithADouble* fSetCustomMatOmassfractCmd;
        G4UIcmdWithADouble* fSetCustomMatNamassfractCmd;
        G4UIcmdWithADouble* fSetCustomMatPmassfractCmd;
        G4UIcmdWithADouble* fSetCustomMatSmassfractCmd;
        G4UIcmdWithADouble* fSetCustomMatClmassfractCmd;
        G4UIcmdWithADouble* fSetCustomMatKmassfractCmd;
        G4UIcmdWithADouble* fSetCustomMatCamassfractCmd;
        
    G4UIcmdWithAnInteger* fPhantomMaterialCmd;        
    G4UIcmdWithADoubleAndUnit* fPhantomDiameterCmd;
    G4UIcmdWithADoubleAndUnit* fPhantomHeightCmd;
    G4UIcmdWithADoubleAndUnit* fPhantomZCmd;
    
    G4UIcmdWithADouble* fSetComp0Cmd;
    G4UIcmdWithADouble* fSetComp1Cmd;
    G4UIcmdWithADouble* fSetComp2Cmd;
    G4UIcmdWithADouble* fSetComp3Cmd;   
    
    G4UIcmdWithADouble* fThetaSetupCmd;   
    
    G4UIcmdWithABool* fSetSlitsCmd;
    G4UIcmdWithADoubleAndUnit* fSlit1ThicknessCmd;
    G4UIcmdWithADoubleAndUnit* fSlit2ThicknessCmd;
    G4UIcmdWithADoubleAndUnit* fSlit3ThicknessCmd;
    G4UIcmdWithADoubleAndUnit* fSlit4ThicknessCmd;
    G4UIcmdWithADoubleAndUnit* fSlit1DistanceCmd;
    G4UIcmdWithADoubleAndUnit* fSlit2DistanceCmd;
    G4UIcmdWithADoubleAndUnit* fSlit3DistanceCmd;
    G4UIcmdWithADoubleAndUnit* fSlit4DistanceCmd;
    G4UIcmdWithADoubleAndUnit* fSlit1xApertureCmd;
    G4UIcmdWithADoubleAndUnit* fSlit2xApertureCmd;
    G4UIcmdWithADoubleAndUnit* fSlit3xApertureCmd;
    G4UIcmdWithADoubleAndUnit* fSlit4xApertureCmd; 
    G4UIcmdWithADoubleAndUnit* fSlit1yApertureCmd;
    G4UIcmdWithADoubleAndUnit* fSlit2yApertureCmd;
    G4UIcmdWithADoubleAndUnit* fSlit3yApertureCmd;
    G4UIcmdWithADoubleAndUnit* fSlit4yApertureCmd; 
    
        G4UIcmdWithADoubleAndUnit* fDetectorThicknessCmd;   
    G4UIcmdWithADoubleAndUnit* fDetectorSizeCmd;
    G4UIcmdWithADoubleAndUnit* fDetectorDistanceCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

