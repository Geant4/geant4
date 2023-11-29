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
/// \file SteppingAction.hh
/// \brief Definition of the SteppingAction class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef SteppingAction_H
#define SteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"
#include <array>

class Run;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SteppingAction : public G4UserSteppingAction {
  public:  
    SteppingAction();
    ~SteppingAction() override = default;
  
    void UserSteppingAction( const G4Step* ) override;
    // This is the main method where the step lengths of particles inside
    // the scoring shell are collected, and then the corresponding fluences
    // are filled up in the Run object where they are stored (and then
    // printed out at the end of the Run).
    // (For simplicity and brevity, we avoid histograms and compute instead
    //  some statistics ourself, which will be print-out at the end of the run.)

    void Initialize();
    // This method is called by RunAction::BeginOfRunAction for the
    // initialization of the stepping-action at the beginning of each Run.
    // This is necessary because different runs can have different primary particle
    // types, kinetic energies, and detector configurations.

    void SetRunPointer( Run* inputValue = nullptr ) { fRunPtr = inputValue; }
    // This method is called by RunAction::BeginOfRunAction for providing to the
    // stepping-action the pointer to the run object at the beginning of each Run.
    // This pointer is then used to pass the information collected by the stepping-action
    // to the run object.

    G4double GetCubicVolumeScoringShell() const { return fCubicVolumeScoringShell; }
    // The cubic-volume of the scoring shell is needed to get the fluence from the
    // sum of step lengths inside that scoring shell.
  
    static const G4int fkNumberKinematicRegions = 3;  // all, below 20 MeV, above 20 MeV
    static const G4int fkNumberScoringPositions = 2;  // forward, backward (hemisphere with
                                                    // respect to the primary particle direction)
    static const G4int fkNumberParticleTypes = 11;    // all, e, gamma, mu, nu, pi, n, p, ions,
                                                    // other-mesons, other-baryons
    static const G4int fkNumberCombinations =
      fkNumberKinematicRegions*fkNumberScoringPositions*fkNumberParticleTypes;
    static const std::array< G4String, fkNumberKinematicRegions > fkArrayKinematicRegionNames;
    static const std::array< G4String, fkNumberScoringPositions > fkArrayScoringPositionNames;
    static const std::array< G4String, fkNumberParticleTypes > fkArrayParticleTypeNames;
    static G4int GetIndex( const G4int iKinematicRegion, const G4int iScoringPosition,
                           const G4int iParticleType );
 
  private:   
    Run* fRunPtr;  // Pointer to the Run object
    G4int fPrimaryParticleId;
    G4double fPrimaryParticleEnergy;
    G4ThreeVector fPrimaryParticleDirection;
    G4String fTargetMaterialName;
    G4bool fIsFirstStepOfTheEvent;
    G4bool fIsFirstStepInTarget;
    G4bool fIsFirstStepInScoringShell;
    G4double fCubicVolumeScoringShell;
  
    std::array< G4double, fkNumberCombinations > fArraySumStepLengths;
    // Array to collect the sum of step lengths in the scoring shell for the whole run,
    // according to the various cases (kinematical region, scoring position and particle type).
    // Note that the fluence in a scoring volume is defined as sum of step lengths
    // in that scoring volume divided by the cubic-volume of that scoring volume.
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
