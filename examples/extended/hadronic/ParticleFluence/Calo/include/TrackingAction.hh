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
/// \file TrackingAction.hh
/// \brief Definition of the TrackingAction class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef TrackingAction_h 
#define TrackingAction_h 1

#include "globals.hh"
#include "G4UserTrackingAction.hh"
#include <array>

class Run;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TrackingAction : public G4UserTrackingAction {
  // We are using this class to monitor the average multiplicity, the average
  // kinetic energy, and the average total energy flow (i.e. sum of the
  // kinetic energies) of different particle types as they are produced
  // inside the calorimeter.
  // The aim is then to try to correlate some changes in these (more primitive)
  // quantities with the observed changes in the (more indirect and complex)
  // particle fluences.
  public:
    TrackingAction();
    ~TrackingAction() override = default;
    void PreUserTrackingAction( const G4Track* ) override;
    void PostUserTrackingAction( const G4Track* ) override;

    void Initialize();
    // This method is called by RunAction::BeginOfRunAction for the
    // initialization of the tracking-action at the beginning of each Run.

    void SetRunPointer( Run* inputValue = nullptr ) { fRunPtr = inputValue; }
    // This method is called by RunAction::BeginOfRunAction for providing to the
    // tracking-action the pointer to the run object at the beginning of each Run.
    // This pointer is then used to pass the information collected by the tracking-action
    // to the run object.
  
    static const G4int fkNumberScoringVolumes = 1;    // calorimeter
    static const G4int fkNumberKinematicRegions = 3;  // all, below 20 MeV, above 20 MeV
    static const G4int fkNumberParticleTypes = 11;    // all, e, gamma, mu, nu, pi, n, p, ions,
                                                    // other-mesons, other-baryons
    static const G4int fkNumberCombinations =
      fkNumberScoringVolumes*fkNumberKinematicRegions*fkNumberParticleTypes;
    static const std::array< G4String, fkNumberScoringVolumes > fkArrayScoringVolumeNames;
    static const std::array< G4String, fkNumberKinematicRegions > fkArrayKinematicRegionNames;
    static const std::array< G4String, fkNumberParticleTypes > fkArrayParticleTypeNames;
    static G4int GetIndex( const G4int iScoringVolume, const G4int iKinematicRegion,
                           const G4int iParticleType );

  private:
    Run* fRunPtr;  // Pointer to the Run object
  
    std::array< G4int,    fkNumberCombinations > fArrayMultiplicities;
    std::array< G4double, fkNumberCombinations > fArraySumKineticEnergies;
    // Keep record of the fkNumber of particles and their kinetic energy at production,
    // according to the particle type and their kinetic energy range (below/above 20 MeV).
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
