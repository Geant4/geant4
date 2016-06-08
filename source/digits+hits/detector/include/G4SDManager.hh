// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SDManager.hh,v 1.2.4.1 1999/12/07 20:47:41 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//

#ifndef G4SDManager_h
#define G4SDManager_h 1

#include "globals.hh"
#include "G4SDStructure.hh"
#include "G4HCtable.hh"
class G4VHitsCollection;
class G4VSensitiveDetector;
class G4HCofThisEvent;
class G4SDmessenger;

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
  public: // with description
      static G4SDManager* GetSDMpointer();
      // Returns the pointer to the singleton object.
  public:
      static G4SDManager* GetSDMpointerIfExist();

  protected:
      G4SDManager();

  public:
      ~G4SDManager();

  public: // with description
      void AddNewDetector(G4VSensitiveDetector*aSD); 
      //  Registors the user's sensitive detector. This method must be invoked
      // when the user construct his/her sensitive detector.
      void Activate(G4String dName, G4bool activeFlag);
      //  Activate/inactivate the registered sensitive detector. For the inactivated
      // detectors, hits collections will not be stored to the G4HCofThisEvent object.
      G4int GetCollectionID(G4String colName);
      G4int GetCollectionID(G4VHitsCollection * aHC);
      //  These two methods return the ID number of the sensitive detector.

  public:
      G4VSensitiveDetector* FindSensitiveDetector(G4String dName);
      G4HCofThisEvent* PrepareNewEvent();
      void TerminateCurrentEvent(G4HCofThisEvent* HCE);

  private: 
      static G4SDManager * fSDManager;
      G4SDStructure * treeTop;
      G4int verboseLevel;
      G4HCtable* HCtable;
      G4SDmessenger* theMessenger;

  public:
      inline void SetVerboseLevel(G4int vl) 
      { 
        verboseLevel = vl; 
        treeTop->SetVerboseLevel(vl);
      }
      inline G4SDStructure* GetTreeTop() const
      { return treeTop; }
      inline void ListTree() const
      { treeTop->ListTree(); }
      inline G4int GetCollectionCapacity() const
      { return HCtable->entries(); }
      inline G4HCtable* GetHCtable() const
      { return HCtable; }

};




#endif

