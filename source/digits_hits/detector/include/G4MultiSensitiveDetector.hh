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
//
#ifndef G4MULTISENSITIVEDETECTOR_H
#define G4MULTISENSITIVEDETECTOR_H

#include "G4VSensitiveDetector.hh"

#include <vector>

// class description:
// This class allows to assign multiple sensitive detectors to a single
// logical-volume.
// SDs are added to this proxy and an instance of the proxy is assigned
// to the logical volume. Calls to SD methods are forwarded to ALL
// user-defined SD that are added.

class G4MultiSensitiveDetector : public G4VSensitiveDetector
{
 public:
  using sds_t = std::vector<G4VSensitiveDetector*>;
  using sdsConstIter = sds_t::const_iterator;

 public:
  using G4VSensitiveDetector::G4VSensitiveDetector;
  ~G4MultiSensitiveDetector() override = default;

  G4MultiSensitiveDetector(const G4MultiSensitiveDetector& rhs) = default;
  G4MultiSensitiveDetector& operator=(const G4MultiSensitiveDetector& rhs) = default;

 public:
  // interface from G4VSensitiveDetector starts here.
  // See G4VSensitiveDetector for documentation.
  // All these methods forward the call to each of the SD
  // attached to this proxy.
  void Initialize(G4HCofThisEvent*) override;
  void EndOfEvent(G4HCofThisEvent*) override;
  void clear() override;
  void DrawAll() override;
  void PrintAll() override;

  // Return clone of this detector
  // Requires all held SDs to be cloneable
  G4VSensitiveDetector* Clone() const override;

  G4VSensitiveDetector* GetSD(const int i) const { return fSensitiveDetectors[i]; }

  sds_t::size_type GetSize() const { return fSensitiveDetectors.size(); }
  sdsConstIter GetBegin() const { return fSensitiveDetectors.begin(); }
  sdsConstIter GetEnd() const { return fSensitiveDetectors.end(); }
  void ClearSDs() { fSensitiveDetectors.clear(); }
  void AddSD(G4VSensitiveDetector* sd) { fSensitiveDetectors.push_back(sd); }

 protected:
  // The return value is an AND of the called SDs return values.
  // This method will call the "Hit(G4Step*)" method of all
  // added SDs. Note that the ROhist of this method is not used
  G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist) override;

  // The following method does not have a meaning for this concrete class
  G4int GetCollectionID(G4int i) final;

 private:
  sds_t fSensitiveDetectors;
};

#endif  // G4MULTISENSITIVEDETECTOR_H
