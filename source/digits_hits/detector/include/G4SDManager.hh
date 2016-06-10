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
// $Id: G4SDManager.hh 81087 2014-05-20 15:44:27Z gcosmo $
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
      G4VSensitiveDetector* FindSensitiveDetector(G4String dName, G4bool warning = true);
      G4HCofThisEvent* PrepareNewEvent();
      void TerminateCurrentEvent(G4HCofThisEvent* HCE);
      void AddNewCollection(G4String SDname,G4String DCname);


  private: 
      static G4ThreadLocal G4SDManager * fSDManager;
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
private:
    //Disable copy constructor and assignment operator
    G4SDManager( const G4SDManager& );
    G4SDManager& operator=(const G4SDManager&);

};




#endif

