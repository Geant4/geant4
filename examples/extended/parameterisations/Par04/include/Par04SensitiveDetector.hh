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
#ifndef PAR04SENSITIVEDETECTOR_HH
#define PAR04SENSITIVEDETECTOR_HH

#include <CLHEP/Units/SystemOfUnits.h>     // for m, pi
#include <G4String.hh>                     // for G4String
#include <G4Types.hh>                      // for G4bool, G4int
#include "G4SystemOfUnits.hh"              // for m
#include "G4ThreeVector.hh"                // for G4ThreeVector
#include "G4VFastSimSensitiveDetector.hh"  // for G4VFastSimSensitiveDetector
#include "G4VSensitiveDetector.hh"         // for G4VSensitiveDetector
#include "Par04Hit.hh"                     // for Par04Hit (ptr only), Par04...
class G4FastHit;
class G4FastTrack;
class G4HCofThisEvent;
class G4Step;
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

class Par04SensitiveDetector
  : public G4VSensitiveDetector
  , public G4VFastSimSensitiveDetector
{
 public:
  Par04SensitiveDetector(G4String aName);
  Par04SensitiveDetector(G4String aName, G4ThreeVector aNbOfCells, G4ThreeVector aNSizeOfCells);
  virtual ~Par04SensitiveDetector();
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
  Par04Hit* RetrieveAndSetupHit(G4ThreeVector aPosition);

 private:
  /// Collection of hits
  Par04HitsCollection* fHitsCollection = nullptr;
  /// ID of collection of hits
  G4int fHitCollectionID = -1;
  /// Number of mesh readout cells in cylindrical coordinates
  G4ThreeVector fMeshNbOfCells = { 10, 10, 10 };
  /// Size of mesh readout cells in cylindrical coordinates.
  G4ThreeVector fMeshSizeOfCells = { 1 * m, 2 * CLHEP::pi / 10., 1 * m };
  /// Retrieved once per event: position of entering particle
  G4ThreeVector fEntrancePosition = { -1, -1, -1 };
  /// Retrieved once per event: direction of entering particle
  G4ThreeVector fEntranceDirection = { -1, -1, -1 };
};

#endif /* PAR04SENSITIVEDETECTOR_HH */
