// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SDManager.hh,v 2.1 1998/07/12 02:53:19 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
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

class G4SDManager 
{
  public:
      static G4SDManager* GetSDMpointer();
      static G4SDManager* GetSDMpointerIfExist();

  protected:
      G4SDManager();

  public:
      ~G4SDManager();
      void AddNewDetector(G4VSensitiveDetector*aSD); 
      G4HCofThisEvent* PrepareNewEvent();
      void TerminateCurrentEvent(G4HCofThisEvent* HCE);
      void Activate(G4String dName, G4bool activeFlag);
      G4int GetCollectionID(G4String colName);
      G4int GetCollectionID(G4VHitsCollection * aHC);
      G4VSensitiveDetector* FindSensitiveDetector(G4String dName);

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

