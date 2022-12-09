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
/// \file Run.hh
/// \brief Definition of the Run class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "G4ThreeVector.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh" 
#include <array>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run {
  // This class accumulates relevant quantities related to particle fluence collected during
  // the run.
  // ( Note: these information are provided via calls of accessor methods of this Run class
  //         made by SteppingAction::UserSteppingAction
  //         and     TrackingAction::PreUserTrackingAction. )
  // At the end of a run, the  PrintInfo  method is called by the run-action to print out
  // some summary information about these quantities.
  // In multithreaded (MT) mode, an object of this class is filled up for each working thread,
  // and then merged (automatically by the Geant4 kernel) into another object (of this class)
  // owned by the master class; the  PrintInfo  method is then called only for the latter run
  // object.
  // Note that, for simplicity and brevity, we avoid histograms and print-out instead some
  // statistics (compute by ourself) at the end of the run.  
  public:
    Run();
    ~Run() override = default;
  
    void RecordEvent( const G4Event* anEvent ) override;
    // This method is called automatically by the Geant4 kernel (not by the user!) at the end
    // of each event. In the case of multithreaded mode, it is called only for the working thread
    // that handled that event.

    void Merge( const G4Run* aRun ) override;
    // This method is called automatically by the Geant4 kernel (not by the user!) only in the
    // case of multithreaded mode and only for working threads.
  
    void PrintInfo() const;
    // This method is called by RunAction::EndOfRunAction : in the case of multithreaded mode,
    // only the master thread calls it.

    void SetPrimaryParticleId( const G4int inputValue ) { fPrimaryParticleId = inputValue; }
    void SetPrimaryParticleEnergy( const G4double inputValue )
      { fPrimaryParticleEnergy = inputValue; }
    void SetPrimaryParticleDirection( const G4ThreeVector &inputValue )
      { fPrimaryParticleDirection = inputValue; }
    void SetAbsorberMaterialName( const G4String &inputValue )
      { fAbsorberMaterialName = inputValue; }
    void SetActiveMaterialName( const G4String &inputValue ) { fActiveMaterialName = inputValue; }
    void SetCubicVolumeScoringUpDown( const G4double inputValue )
      { fCubicVolumeScoringUpDown = inputValue; }
    void SetCubicVolumeScoringSide( const G4double inputValue )
      { fCubicVolumeScoringSide = inputValue; }
    G4int GetPrimaryParticleId() const { return fPrimaryParticleId; }
    G4double GetPrimaryParticleEnergy() const { return fPrimaryParticleEnergy; }
    G4ThreeVector GetPrimaryParticleDirection() const { return fPrimaryParticleDirection; }
    G4String GetAbsorberMaterialName() const { return fAbsorberMaterialName; }
    G4String GetActiveMaterialName() const { return fActiveMaterialName; }
    G4double GetCubicVolumeScoringUpDown() const { return fCubicVolumeScoringUpDown; }
    G4double GetCubicVolumeScoringSide() const { return fCubicVolumeScoringSide; }

    void SetSteppingArray( const std::array< G4double,
                           SteppingAction::fkNumberCombinations >& inputArray );
    std::array< G4double, SteppingAction::fkNumberCombinations > GetSteppingArray() const
    { return fSteppingArray; }
    // Accessor methods useful to transfer information collected by the stepping-action
    // into this Run class

    void SetTrackingArray1( const std::array< G4int,
                            TrackingAction::fkNumberCombinations >& inputArray );
    std::array< G4int, TrackingAction::fkNumberCombinations > GetTrackingArray1() const
    { return fTrackingArray1; }
    void SetTrackingArray2( const std::array< G4double,
                            TrackingAction::fkNumberCombinations >& inputArray );
    std::array< G4double, TrackingAction::fkNumberCombinations > GetTrackingArray2() const
    { return fTrackingArray2; }
    // Accessor methods useful to transfer information collected by the tracking-action
    // into this Run class

  private:  
    G4int fNumEvents;
    G4int fPrimaryParticleId;
    G4double fPrimaryParticleEnergy;
    G4ThreeVector fPrimaryParticleDirection;
    G4String fAbsorberMaterialName;
    G4String fActiveMaterialName;
    G4double fCubicVolumeScoringUpDown;
    G4double fCubicVolumeScoringSide;
    std::array< G4double, SteppingAction::fkNumberCombinations > fSteppingArray;
    std::array< G4int,    TrackingAction::fkNumberCombinations > fTrackingArray1;
    std::array< G4double, TrackingAction::fkNumberCombinations > fTrackingArray2;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
