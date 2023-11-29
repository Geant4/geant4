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
// G4Step
//
// Class description:
//
// This class represents the Step of a particle being tracked.
// It includes information of:
//     1) List of Step points which compose the Step,
//     2) static information of particle which generated the Step,
//     3) trackID and parent particle ID of the Step,
//     4) termination condition of the Step.

// Authors:
//   Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//   Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
// Revisions:
//   Hisaya Kurashige, 1998-2007
// --------------------------------------------------------------------
#ifndef G4Step_hh
#define G4Step_hh 1

#include <cstdlib>              // Include from 'system'
#include <cmath>                 // Include from 'system'
#include "G4ios.hh"              // Include from 'system'
#include <iomanip>               // Include from 'system'
#include "globals.hh"            // Include from 'global'
#include "G4ThreeVector.hh"      // Include from 'global'
#include "G4VPhysicalVolume.hh"  // Include from 'geometry'
#include "G4StepPoint.hh"        // Include from 'track'
#include "G4StepStatus.hh"       // Include from 'track'
#include "G4TrackVector.hh"      // Include from 'tracking'
#include "G4Profiler.hh"         // Include from 'global'

class G4Polyline;                // Forward declaration.
class G4Track;                   // Forward declaration.

class G4Step
{
  public:
   // the profiler aliases are only used when compiled with GEANT4_USE_TIMEMORY
   using ProfilerConfig = G4ProfilerConfig<G4ProfileType::Step>;

    G4Step();
   ~G4Step();
      // Constructor/Destructor

    G4Step(const G4Step&);
    G4Step& operator=(const G4Step&);
      // Copy Constructor and assignment operator

    G4Track* GetTrack() const;
    void SetTrack(G4Track* value);
      // Current track

    G4StepPoint* GetPreStepPoint() const;
    void SetPreStepPoint(G4StepPoint* value);
    G4StepPoint* ResetPreStepPoint(G4StepPoint* value=nullptr);
      // Pre-Step points
      // If Set method is invoked, the previous StepPoint object is deleted.
      // If Reset method is invoked, the previous StepPoint object is not deleted
      // but its pointer is returned. Thus it's the caller's responsibility to 
      // properly delete it.

    G4StepPoint* GetPostStepPoint() const;
    void SetPostStepPoint(G4StepPoint* value);
    G4StepPoint* ResetPostStepPoint(G4StepPoint* value=nullptr);
      // Post-Step points
      // If Set method is invoked, the previous StepPoint object is deleted.
      // If Reset method is invoked, the previous StepPoint object is not deleted
      // but its pointer is returned. Thus it's the caller's responsibility to 
      // properly delete it.

    G4double GetStepLength() const;
    void SetStepLength(G4double value);
      // Before the end of the AlongStepDoIt loop, StepLength keeps
      // the initial value which is determined by the shortest geometrical Step
      // proposed by a physics process. After finishing the AlongStepDoIt,
      // it will be set equal to 'StepLength' in G4Step

    G4double GetTotalEnergyDeposit() const;
    void SetTotalEnergyDeposit(G4double value);
      // Total energy deposit

    G4double GetNonIonizingEnergyDeposit() const;
    void SetNonIonizingEnergyDeposit(G4double value);
      // Total non-ionizing energy deposit

    G4SteppingControl GetControlFlag() const;
    void SetControlFlag(G4SteppingControl StepControlFlag);
      // Control flag for stepping

    void AddTotalEnergyDeposit(G4double value);
    void ResetTotalEnergyDeposit();
      // Manipulation of total energy deposit

    void AddNonIonizingEnergyDeposit(G4double value);
    void ResetNonIonizingEnergyDeposit();
      // Manipulation of non-ionizing energy deposit

    G4bool IsFirstStepInVolume() const;
    G4bool IsLastStepInVolume() const;

    void SetFirstStepFlag();
    void ClearFirstStepFlag();
    void SetLastStepFlag();
    void ClearLastStepFlag();
      // Get/Set/Clear flag for initial/last step
      // NOTE: flags are not used

    G4ThreeVector GetDeltaPosition() const;
    G4double GetDeltaTime() const;
      // Difference of position, time, momentum and energy

    G4ThreeVector GetDeltaMomentum() const;
    G4double GetDeltaEnergy() const;
      // These methods will be deleted
      // NOTE: use  GetTotalEnergyDeposit() to obtain energy loss in material

    void InitializeStep(G4Track* aValue);
      // Initialize contents of G4Step

    void UpdateTrack();
      // Update track by using G4Step information

    void CopyPostToPreStepPoint();
      // Copy PostStepPoint to PreStepPoint

    G4Polyline* CreatePolyline() const;
      // For visualization

    inline void SetPointerToVectorOfAuxiliaryPoints(std::vector<G4ThreeVector>* vec);
    inline std::vector<G4ThreeVector>* GetPointerToVectorOfAuxiliaryPoints() const;
      // Auxiliary points modifiers

  // --- Secondary buckets ---

    std::size_t GetNumberOfSecondariesInCurrentStep() const;
      // Secondaries in the current step

    const std::vector<const G4Track*>* GetSecondaryInCurrentStep() const;

    const G4TrackVector* GetSecondary() const;
    G4TrackVector* GetfSecondary();
    G4TrackVector* NewSecondaryVector();
      // NOTE: Secondary bucket of the Step contains
      //       all secondaries during tracking the current track
      //       (i.e. NOT secondaries produced in the current step)
      // all these methods give same object (i.e. G4TrackVector  )
      // but 2nd one will create bucket in addition

    void DeleteSecondaryVector();
      // Just delete secondary bucket
      // NOTE: G4Track objects inside the bucket are not deleted

    void SetSecondary(G4TrackVector* value);
      // Add secondary tracks to the bucket

  protected:

    G4double fTotalEnergyDeposit = 0.0;
      // Accumulated total energy deposit in the current Step

    G4double fNonIonizingEnergyDeposit = 0.0;
     // Accumulated non-ionizing energy deposit in the current Step

  private:

    G4StepPoint* fpPreStepPoint = nullptr;
    G4StepPoint* fpPostStepPoint = nullptr;
    G4double fStepLength = 0.0;
      // Step length which may be updated at each invocation of
      // AlongStepDoIt and PostStepDoIt

    G4Track* fpTrack = nullptr;

    G4SteppingControl fpSteppingControlFlag = NormalCondition;
      // A flag to control SteppingManager behavior from process

    G4bool fFirstStepInVolume = false;
    G4bool fLastStepInVolume = false;
      // Flags for initial/last step

    G4TrackVector* fSecondary = nullptr;
      // Secondary bucket implemented by using  std::vector of G4Track*

    std::size_t nSecondaryByLastStep = 0;
      // number of secondaries which have been created by the last step

    std::vector<const G4Track*>* secondaryInCurrentStep = nullptr;

    std::vector<G4ThreeVector>* fpVectorOfAuxiliaryPointsPointer = nullptr;
};

#include "G4Step.icc"

#endif
