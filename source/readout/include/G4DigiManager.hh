// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DigiManager.hh,v 1.1 1999-01-07 16:14:13 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4DigiManager_h
#define G4DigiManager_h 1

#include "globals.hh"
class G4Event;
#include "G4VDigitizerModule.hh"
class G4VHitsCollection;
class G4VDigiCollection;
class G4DMmessenger;
#include "G4DCtable.hh"
class G4RunManager;
class G4SDManager;
#include <rw/tpordvec.h>

class G4DigiManager 
{
  public:
      static G4DigiManager* GetDMpointer();
      static G4DigiManager* GetDMpointerIfExist();

  protected:
      G4DigiManager();

  public:
      ~G4DigiManager();
      void AddNewModule(G4VDigitizerModule* DM); 
      void Digitize(G4String mName);
      G4VDigitizerModule* FindDigitizerModule(G4String mName);
      const G4VHitsCollection* GetHitsCollection(G4int HCID, G4int eventID = 0);
      const G4VDigiCollection* GetDigiCollection(G4int DCID, G4int eventID = 0);
      G4int GetHitsCollectionID(G4String HCname);
      G4int GetDigiCollectionID(G4String DCname);
      void SetDigiCollection(G4int DCID, G4VDigiCollection* aDC);
      void SetVerboseLevel(G4int vl);
      void List() const;

  private: 
      static G4DigiManager * fDManager;
      G4int verboseLevel;
      RWTPtrOrderedVector<G4VDigitizerModule> DMtable;
      G4DCtable* DCtable;
      G4DMmessenger* theMessenger;
      G4RunManager* runManager;
      G4SDManager* SDManager;

  public:
      inline G4int GetVerboseLevel() const
      { return verboseLevel; }
      inline G4int GetCollectionCapacity() const
      { return DCtable->entries(); }
      inline G4int GetModuleCapacity() const
      { return DMtable.entries(); }
      inline G4DCtable* GetDCtable() const
      { return DCtable; }
      inline void RestoreDCtable(G4DCtable* dc)
      { 
        if(DCtable) delete DCtable;
        DCtable = dc;
      }

};




#endif

