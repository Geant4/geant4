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
//

#ifndef G4VPrimitiveScorer_h
#define G4VPrimitiveScorer_h 1

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;
#include "G4MultiFunctionalDetector.hh"
#include "G4VSDFilter.hh"
#include "globals.hh"

// class description:
//
//  This is the base class of the sensitive detector which owns
// only one hits collection.
//  A concrete class object derived from this base class can be
// used either as a sensitive detector or to be registered to
// G4MultiFunctionalDetector to define multiple functionalities.
//
//

class G4VPrimitiveScorer
{
  friend class G4MultiFunctionalDetector;

 public:
  G4VPrimitiveScorer(G4String name, G4int depth = 0);
  virtual ~G4VPrimitiveScorer() = default;

  // This method returns the ID of its hitsCollection. This mehod
  // gives valid value only after it is registered to G4MultiFunctionalDetector
  // and the G4MultiFunctionalDetector is registered to G4SDManager.
  G4int GetCollectionID(G4int);

  // These five methods are exactly identical to those in G4VSensitiveDetector.
  // These methods are invoked by G4SDManager through G4MultiFunctionalDetector.
  virtual void Initialize(G4HCofThisEvent*);
  virtual void EndOfEvent(G4HCofThisEvent*);
  virtual void clear();
  virtual void DrawAll();
  virtual void PrintAll();

  void SetUnit(const G4String& unit) { unitName = unit; }
  const G4String& GetUnit() const { return unitName; }
  G4double GetUnitValue() const { return unitValue; }

  // Set/Get methods
  inline void SetMultiFunctionalDetector(G4MultiFunctionalDetector* d) { detector = d; }
  inline G4MultiFunctionalDetector* GetMultiFunctionalDetector() const { return detector; }
  inline G4String GetName() const { return primitiveName; }
  inline void SetFilter(G4VSDFilter* f) { filter = f; }
  inline G4VSDFilter* GetFilter() const { return filter; }
  inline void SetVerboseLevel(G4int vl) { verboseLevel = vl; }
  inline G4int GetVerboseLevel() const { return verboseLevel; }

  inline void SetNijk(G4int i, G4int j, G4int k)
  {
    fNi = i;
    fNj = j;
    fNk = k;
  }

 protected:
  // Get the solid at current depth, ensuring it's correct by
  //   calling a parameterisation is called if it's that volume type
  G4VSolid* ComputeSolid(G4Step* aStep, G4int replicaIdx);

  // Same as above -- using stored replica number
  G4VSolid* ComputeCurrentSolid(G4Step* aStep);

  // This is the method must be implemented in each concrete class.
  virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*) = 0;

  // This is a function mapping from copy number(s) to an index of
  // the hit collection. In the default implementation, just the
  // copy number of the physical volume is taken.
  virtual G4int GetIndex(G4Step*);

  void CheckAndSetUnit(const G4String& unit, const G4String& category);

 protected:
  G4String primitiveName;
  G4MultiFunctionalDetector* detector{nullptr};
  G4VSDFilter* filter{nullptr};
  G4int verboseLevel{0};
  G4int indexDepth;
  G4String unitName{"NoUnit"};
  G4double unitValue{1.0};
  G4int fNi{0}, fNj{0}, fNk{0};  // used for 3D scorers

 private:
  inline G4bool HitPrimitive(G4Step* aStep, G4TouchableHistory* ROhis)
  {
    if (filter != nullptr) {
      if (! (filter->Accept(aStep))) return false;
    }
    return ProcessHits(aStep, ROhis);
  }
};

#endif
