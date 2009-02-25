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
#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SteppingAction;
class DetectorConstruction;
class PrimaryGeneratorAction;

class G4Run;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
public:

  RunAction(DetectorConstruction*, PrimaryGeneratorAction*);
  virtual ~RunAction();

  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);

  void SumEnergy1(G4int k, G4double de){energyDeposit1[k] += de; };  	
  void SumEnergy2(G4int k, G4double de){energyDeposit2[k] += de; };  	
  void SumEnergy3(G4int k, G4double de){energyDeposit3[k] += de; };  	

  void AddEnergy1(G4double de) {energyDepositRun1 += de; n_steps++;};
  void AddEnergy2(G4double de) {energyDepositRun2 += de; n_steps++;};
  void AddEnergy3(G4double de) {energyDepositRun3 += de; n_steps++;};
    
private:
  DetectorConstruction*   detector;
  PrimaryGeneratorAction* primary;

  G4int n_steps;

  G4String matName1     ;
  G4String matName2     ;
  G4String matName3     ;

  G4double              energyDeposit1[110];
  G4double              energyDeposit2[110];
  G4double              energyDeposit3[110];
  G4double              normalizedvalue1[110];
  G4double              normalizedvalue2[110];
  G4double              normalizedvalue3[110];
  G4double MFP1;        
  G4double MFP2;        
  G4double MFP3;        
  G4double density1;        
  G4double density2;        
  G4double density3;        
  G4double energyDepositRun1; 
  G4double energyDepositRun2;
  G4double energyDepositRun3;
  G4String asciiFileName;
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

