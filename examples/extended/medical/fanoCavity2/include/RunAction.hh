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
// $Id: RunAction.hh,v 1.2 2007-10-31 16:16:20 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "ProcessesCount.hh"
#include "globals.hh"

class G4Run;
class G4Material;
class DetectorConstruction;
class PrimaryGeneratorAction;
class HistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:
    RunAction(DetectorConstruction*, PrimaryGeneratorAction*,
              HistoManager*);
   ~RunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);

    void CountProcesses(G4String);
    
    void SurveyConvergence(G4int);
    
    void FlowInCavity(G4int k, G4double e) { EnerFlowCavity[k] += e;  
                                             PartFlowCavity[k]++;};
					     
    void AddEdepCavity(G4double de) { EdepCavity += de; EdepCavity2 += de*de;
                                      nbEventCavity++;};
    void AddTrakCavity(G4double dt) { trkSegmCavity += dt;};
        
    void StepInWall   (G4double s)  { stepWall += s; stepWall2 += s*s; 
                                      nbStepWall++;};
    void StepInCavity (G4double s)  { stepCavity += s; stepCavity2 += s*s; 
                                      nbStepCavity++;};
    
  private:
    DetectorConstruction*   detector;
    PrimaryGeneratorAction* kinematic;
    ProcessesCount*         ProcCounter;
    HistoManager*           histoManager;
    
    G4long                  PartFlowCavity[2];
    G4double                EnerFlowCavity[2];        
    G4double                EdepCavity, EdepCavity2;
    G4double                trkSegmCavity;
    G4long                  nbEventCavity;
                
    G4double                stepWall,   stepWall2;
    G4double                stepCavity, stepCavity2;    
    G4long                  nbStepWall, nbStepCavity;
    
  private:
    G4double                energyGun;  
    G4double                massWall;
    G4double                massCavity;  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

