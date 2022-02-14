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
#ifndef PAR03SENSITIVEDETECTOR_HH
#define PAR03SENSITIVEDETECTOR_HH

#include "Par03Hit.hh"

#include "G4VFastSimSensitiveDetector.hh"
#include "G4VSensitiveDetector.hh"

class G4HCofThisEvent;
class G4TouchableHistory;

/**
 * @brief Sensitive detector.
 *
 * Describes how to store the energy deposited within the detector.
 * It derives from two classes: G4VSensitiveDetector and
 * G4VFastSimSensitiveDetector. Addition of G4VFastSimSensitiveDetector is
 * necessary in order to handle the energy deposits from the fast simulation.
 *
 * Two ProcessHits() methods are introduced to handle energy deposited from full
 * (detailed) simulation, and from fast simulation. The common part is handled
 * by RetrieveAdnSetupHit() method.
 *
 */

class Par03SensitiveDetector
  : public G4VSensitiveDetector
  , public G4VFastSimSensitiveDetector
{
 public:
  Par03SensitiveDetector(G4String aName);
  Par03SensitiveDetector(G4String aName, G4int aNumLayers, G4int aNumPhi,
                         G4int aNumRho);
  virtual ~Par03SensitiveDetector();
  /// Create hit collection
  virtual void Initialize(G4HCofThisEvent* HCE) final;
  /// Process energy deposit from the full simulation.
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* aROhist) final;
  /// Process energy deposit from the fast simulation.
  virtual G4bool ProcessHits(const G4FastHit* aHit, const G4FastTrack* aTrack,
                             G4TouchableHistory* aROhist) final;
  /// Process energy deposit - common part for full and fast simulation
  /// It is invoked from ProcessHits() methods, and sets basic hit properties
  /// (position, etc.), common for hit from fast and full simulation.
  Par03Hit* RetrieveAndSetupHit(G4TouchableHistory* aTouchable);

 private:
  /// Collection of hits
  Par03HitsCollection* fHitsCollection = nullptr;
  /// ID of collection of hits
  G4int fHitCollectionID = -1;
  /// Number of readout cells along z axis
  G4int fCellNoZ = 10;
  /// Number of readout cells along radius of cylinder
  G4int fCellNoRho = 10;
  /// Number of readout cells along azimuthal angle
  G4int fCellNoPhi = 10;
};

#endif /* PAR03SENSITIVEDETECTOR_HH */
