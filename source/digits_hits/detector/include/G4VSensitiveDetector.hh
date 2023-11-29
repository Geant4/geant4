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

#ifndef G4VSensitiveDetector_h
#define G4VSensitiveDetector_h 1

#include "G4CollectionNameVector.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4VHit.hh"
#include "G4VReadOutGeometry.hh"
#include "G4VSDFilter.hh"

// class description:
//
//  This is the abstract base class of the sensitive detector. The user's
// sensitive detector which generates hits must be derived from this
// class.
//  In the derived class constructor, name(s) of hits collection(s) which
// are made by the sensitive detector must be set to "collectionName" string
// vector.

class G4VSensitiveDetector
{
 public:
  // Constructors. The user's concrete class must use one of these constructors
  // by the constructor initializer of the derived class. The name of
  // the sensitive detector must be unique.
  explicit G4VSensitiveDetector(G4String name);
  G4VSensitiveDetector(const G4VSensitiveDetector& right);
  G4VSensitiveDetector& operator=(const G4VSensitiveDetector& right);
  virtual ~G4VSensitiveDetector() = default;

  G4bool operator==(const G4VSensitiveDetector& right) const;
  G4bool operator!=(const G4VSensitiveDetector& right) const;

  //  These two methods are invoked at the begining and at the end of each
  // event. The hits collection(s) created by this sensitive detector must
  // be set to the G4HCofThisEvent object at one of these two methods.
  virtual void Initialize(G4HCofThisEvent*) {}
  virtual void EndOfEvent(G4HCofThisEvent*) {}

  //  This method is invoked if the event abortion is occured. Hits collections
  // created but not beibg set to G4HCofThisEvent at the event should be
  // deleted. Collection(s) which have already set to G4HCofThisEvent will be
  // deleted automatically.
  virtual void clear() {}

  virtual void DrawAll() {}
  virtual void PrintAll() {}

  //  This is the public method invoked by G4SteppingManager for generating
  // hit(s). The actual user's implementation for generating hit(s) must be
  // implemented in GenerateHits() virtual protected method. This method
  // MUST NOT be overriden.
  inline G4bool Hit(G4Step* aStep)
  {
    G4TouchableHistory* ROhis = nullptr;
    if (! isActive()) return false;
    if (filter != nullptr) {
      if (! (filter->Accept(aStep))) return false;
    }
    if (ROgeometry != nullptr) {
      if (! (ROgeometry->CheckROVolume(aStep, ROhis))) return false;
    }
    return ProcessHits(aStep, ROhis);
  }

  //  Register the Readout geometry.
  inline void SetROgeometry(G4VReadOutGeometry* value) { ROgeometry = value; }

  //  Register a filter
  inline void SetFilter(G4VSDFilter* value) { filter = value; }

  inline G4int GetNumberOfCollections() const { return G4int(collectionName.size()); }
  inline G4String GetCollectionName(G4int id) const { return collectionName[id]; }
  inline void SetVerboseLevel(G4int vl) { verboseLevel = vl; }
  inline void Activate(G4bool activeFlag) { active = activeFlag; }
  inline G4bool isActive() const { return active; }
  inline G4String GetName() const { return SensitiveDetectorName; }
  inline G4String GetPathName() const { return thePathName; }
  inline G4String GetFullPathName() const { return fullPathName; }
  inline G4VReadOutGeometry* GetROgeometry() const { return ROgeometry; }
  inline G4VSDFilter* GetFilter() const { return filter; }

  virtual G4VSensitiveDetector* Clone() const;

 protected:  // with description
  //  The user MUST implement this method for generating hit(s) using the
  // information of G4Step object. Note that the volume and the position
  // information is kept in PreStepPoint of G4Step.
  //  Be aware that this method is a protected method and it sill be invoked
  // by Hit() method of Base class after Readout geometry associated to the
  // sensitive detector is handled.
  //  "ROhist" will be given only is a Readout geometry is defined to this
  // sensitive detector. The G4TouchableHistory object of the tracking geometry
  // is stored in the PreStepPoint object of G4Step.
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist) = 0;

  //  This is a utility method which returns the hits collection ID of the
  // "i"-th collection. "i" is the order (starting with zero) of the collection
  // whose name is stored to the collectionName protected vector.
  virtual G4int GetCollectionID(G4int i);

 protected:
  //  This protected name vector must be filled at the constructor of the user's
  // concrete class for registering the name(s) of hits collection(s) being
  // created by this particular sensitive detector.
  G4CollectionNameVector collectionName;

  G4String SensitiveDetectorName;  // detector name
  G4String thePathName;  // directory path
  G4String fullPathName;  // path + detector name
  G4int verboseLevel{0};
  G4bool active{true};
  G4VReadOutGeometry* ROgeometry{nullptr};
  G4VSDFilter* filter{nullptr};
};

#endif
