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
// gpaterno, October 2025
//
/// \file PrimaryGeneratorActionMessenger.hh
/// \brief Description of the PrimaryGeneratorActionMessenger class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGeneratorActionMessenger_h
#define PrimaryGeneratorActionMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGeneratorAction;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// messenger for PrimaryGenerator class.

class PrimaryGeneratorActionMessenger: public G4UImessenger
{
public:
    PrimaryGeneratorActionMessenger(PrimaryGeneratorAction*);
    ~PrimaryGeneratorActionMessenger() override;
    
    void SetNewValue(G4UIcommand*, G4String) override;
    
private:   
    PrimaryGeneratorAction* fPrimaryGeneratorAction{nullptr};
    
    G4UIcmdWithABool* fUseGPSCmd{nullptr};
    
    G4UIcmdWithAString* fPrimaryTypeCmd{nullptr}; 
    G4UIcmdWithADoubleAndUnit* fPrimaryEnergyCmd{nullptr}; 
    G4UIcmdWithADouble* fPrimaryRelSigmaEnergyCmd{nullptr}; 
    G4UIcmdWithADoubleAndUnit* fPrimaryXCmd{nullptr};  
    G4UIcmdWithADoubleAndUnit* fPrimaryYCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fPrimaryZCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fPrimaryTCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fPrimaryXpCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fPrimaryYpCmd{nullptr};  
    G4UIcmdWithADoubleAndUnit* fPrimarySigmaXCmd{nullptr}; 
    G4UIcmdWithADoubleAndUnit* fPrimarySigmaYCmd{nullptr}; 
    G4UIcmdWithADoubleAndUnit* fPrimarySigmaZCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fPrimarySigmaTCmd{nullptr}; 
    G4UIcmdWithADoubleAndUnit* fPrimarySigmaXpCmd{nullptr}; 
    G4UIcmdWithADoubleAndUnit* fPrimarySigmaYpCmd{nullptr};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

