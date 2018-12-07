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
/// \file medical/fanoCavity/include/RunAction.hh
/// \brief Definition of the RunAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "ProcessesCount.hh"
#include "globals.hh"

class G4Run;
class G4Material;
class G4ParticleDefinition;
class DetectorConstruction;
class PrimaryGeneratorAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
  public:
    Run(DetectorConstruction*,PrimaryGeneratorAction*);
   ~Run();

  public:
    virtual void Merge(const G4Run*);
    void SurveyConvergence(G4int);

    void EndOfRun();

    void CountProcesses(G4String);
    void SumEsecond(G4double e)    { fEsecondary += e; fEsecondary2 += e*e;
                    fNbSec++;};

    void FlowInCavity(G4int k, G4double e) { fEnerFlowCavity[k] += e;
                        fPartFlowCavity[k]++;};

    void AddEdepCavity(G4double de) { fEdepCavity += de; fEdepCavity2 += de*de;
                                     fNbEventCavity++;};
    void AddTrakCavity(G4double dt) { fTrkSegmCavity += dt;};

    void StepInWall   (G4double s)  { fStepWall += s; fStepWall2 += s*s;
                                        fNbStepWall++;};
    void StepInCavity (G4double s)  { fStepCavity += s; fStepCavity2 += s*s;
                      fNbStepCavity++;};

  private:
    DetectorConstruction*   fDetector;
    PrimaryGeneratorAction* fKinematic;
    ProcessesCount*         fProcCounter;

    G4double                fEsecondary, fEsecondary2;
    G4long                  fNbSec;

    G4long                  fPartFlowCavity[2];
    G4double                fEnerFlowCavity[2];
    G4double                fEdepCavity, fEdepCavity2;
    G4double                fTrkSegmCavity;
    G4long                  fNbEventCavity;

    G4double                fOldEmean, fOldDose;

    G4double                fStepWall,   fStepWall2;
    G4double                fStepCavity, fStepCavity2;
    G4long                  fNbStepWall, fNbStepCavity;

  private:
    G4double                fWallThickness;
    G4double                fWallRadius;
    G4Material*             fMateWall;
    G4double                fDensityWall;

    G4double                fCavityThickness;
    G4double                fCavityRadius;
    G4double                fSurfaceCavity;
    G4double                fVolumeCavity;
    G4Material*             fMateCavity;
    G4double                fDensityCavity;
    G4double                fMassCavity;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

