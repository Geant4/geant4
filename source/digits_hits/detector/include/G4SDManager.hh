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

#ifndef G4SDManager_h
#define G4SDManager_h 1

#include "G4HCtable.hh"
#include "G4SDStructure.hh"
#include "globals.hh"
class G4VHitsCollection;
class G4VSensitiveDetector;
class G4HCofThisEvent;
class G4SDmessenger;

#include "G4VSDFilter.hh"

#include <vector>

// class description:
//
//  This is a singleton class which manages the sensitive detectors.
// The user cannot access to the constructor. The pointer of the
// only existing object can be got via G4SDManager::GetSDMpointer()
// static method. The first invokation of this static method makes
// the singleton object.
//

class G4SDManager
{
 public:
  // Returns the pointer to the singleton object, creating it if not null
  static G4SDManager* GetSDMpointer();

  // Returns current pointer to the singleton object
  // Caller is responsible for checking value against `nullptr`
  static G4SDManager* GetSDMpointerIfExist();

  G4SDManager(const G4SDManager&) = delete;
  G4SDManager& operator=(const G4SDManager&) = delete;
  ~G4SDManager();

  // Register sensitive detector instance
  // This method must be invoked when the user constructs their sensitive detector
  void AddNewDetector(G4VSensitiveDetector* aSD);

  // Activate/inactivate the registered sensitive detector.
  // For the inactivated detectors, hits collections will not be stored to the G4HCofThisEvent
  // object.
  void Activate(G4String dName, G4bool activeFlag);

  // Return ID number of sensitive detector with given name
  G4int GetCollectionID(G4String colName);

  // Return ID number of sensitive detector creating given hits collection
  G4int GetCollectionID(G4VHitsCollection* aHC);

  G4VSensitiveDetector* FindSensitiveDetector(G4String dName, G4bool warning = true);
  G4HCofThisEvent* PrepareNewEvent();
  void TerminateCurrentEvent(G4HCofThisEvent* HCE);
  void AddNewCollection(G4String SDname, G4String DCname);

  inline void SetVerboseLevel(G4int vl)
  {
    verboseLevel = vl;
    treeTop->SetVerboseLevel(vl);
  }
  inline G4SDStructure* GetTreeTop() const { return treeTop; }
  inline void ListTree() const { treeTop->ListTree(); }
  inline G4int GetCollectionCapacity() const { return HCtable->entries(); }
  inline G4HCtable* GetHCtable() const { return HCtable; }

  void RegisterSDFilter(G4VSDFilter* filter);
  void DeRegisterSDFilter(G4VSDFilter* filter);

 protected:
  G4SDManager();

 private:
  void DestroyFilters();

 private:
  static G4ThreadLocal G4SDManager* fSDManager;
  G4SDStructure* treeTop;
  G4int verboseLevel{0};
  G4HCtable* HCtable;
  G4SDmessenger* theMessenger;
  std::vector<G4VSDFilter*> FilterList;
};

#endif
