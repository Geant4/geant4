// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DigiManager.hh,v 1.4 1999-12-15 14:53:51 gunter Exp $
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
#include "g4rw/tpordvec.h"

// class description:
//
//  This is a singleton class which manages the digitizer modules.
// The user cannot access to the constructor. The pointer of the
// only existing object can be got via G4DigiManager::GetDMpointer()
// static method. The first invokation in the program makes the 
// singleton object.
//

class G4DigiManager 
{
  public: // with description
      static G4DigiManager* GetDMpointer();
      // Returns the pointer to the singleton object
  public:
      static G4DigiManager* GetDMpointerIfExist();

  protected:
      G4DigiManager();

  public:
      ~G4DigiManager();

  public: // with description
      void AddNewModule(G4VDigitizerModule* DM); 
      //  Registers the user's digitizer mudule. This method must be invoked when
      // the user construct his/her digitizer module(s).
      void Digitize(G4String mName);
      //  Invokes Digitize() method of specified digitizer module. This is a kind
      // of service method. The user can invoke Digitize() method of a particular
      // module without knowing the pointer of the module object. The argument
      // "mName" is the name of the module, which is defined at the constructor
      // of the concrete digitizer module.
      G4VDigitizerModule* FindDigitizerModule(G4String mName);
      //  Returns the pointer to the digitizer module object with the given name.
      // Null will be returned if the name is not defined.
      const G4VHitsCollection* GetHitsCollection(G4int HCID, G4int eventID = 0);
      const G4VDigiCollection* GetDigiCollection(G4int DCID, G4int eventID = 0);
      //  These two methods return the pointer to the hits and digi collection
      // object, respectively. "HCID" and "DCID" are the ID numbers of hits and 
      // digi collections, which can be obtained vir the next two methods.
      //  If "eventID" is greater than zero, corresponding hits or digi collection
      // of "eventID" prevuois event is returned so that event overlap can be 
      // handled. To do this, necessary number of events must be set to G4RunManager
      // by G4RunManager::SetNumberOfEventsToBeStored() method previously to the
      // event loop.
      G4int GetHitsCollectionID(G4String HCname);
      G4int GetDigiCollectionID(G4String DCname);
      //  Returns the ID number of hits and digi collections, respectively. "HCname"
      // and "DCname" can be the collection name if it is unique, or can be detector
      // or module name and the collection name connected by "/".
      void SetDigiCollection(G4int DCID, G4VDigiCollection* aDC);
      //  This method must exclusively used by the base class of G4VDigitizerModule.
      // To set digi collection, the user must use SetDigiCollection of G4VDigitizerModule.

  public:
      void SetVerboseLevel(G4int vl);
      void List() const;

  private: 
      static G4DigiManager * fDManager;
      G4int verboseLevel;
      G4RWTPtrOrderedVector<G4VDigitizerModule> DMtable;
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

