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
#ifndef PAR04PARALLELFULLSENSITIVEDETECTOR_HH
#define PAR04PARALLELFULLSENSITIVEDETECTOR_HH

#include <CLHEP/Units/SystemOfUnits.h>     // for m, pi
#include <G4String.hh>                     // for G4String
#include <G4Types.hh>                      // for G4bool, G4int
#include "G4SystemOfUnits.hh"              // for m
#include "G4ThreeVector.hh"                // for G4ThreeVector
#include "G4VFastSimSensitiveDetector.hh"  // for G4VFastSimSensitiveDetector
#include "G4VSensitiveDetector.hh"         // for G4VSensitiveDetector
#include "Par04Hit.hh"                     // for Par04Hit (ptr only), Par04...
#include <unordered_map>
class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

/**
 * @brief Sensitive detector.
 *
 *
 */

class Par04ParallelFullSensitiveDetector
  : public G4VSensitiveDetector
{
 public:
  Par04ParallelFullSensitiveDetector(G4String aName);
  Par04ParallelFullSensitiveDetector(G4String aName,
                                     G4int aNbOfLayers, G4int aNbOfSlices, G4int aNbOfRows);
  virtual ~Par04ParallelFullSensitiveDetector();
  /// Create hit collection
  virtual void Initialize(G4HCofThisEvent* HCE) final;
  /// Process energy deposit from the full simulation.
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* aROhist) final;
  virtual void EndOfEvent(G4HCofThisEvent* aHC) final;

 private:
  /// Collection of hits
  Par04HitsCollection* fHitsCollection = nullptr;
  std::unordered_map<G4int, std::unique_ptr<Par04Hit>> fHitsMap;
  /// ID of collection of hits
  G4int fHitCollectionID = -1;
  /// Number of readout cells
  G4int fNbOfLayers = 1;
  G4int fNbOfSlices = 1;
  G4int fNbOfRows = 1;
};

#endif /* PAR04PARALLELSENSITIVEDETECTOR_HH */
