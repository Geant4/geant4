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
// G4PersistencyManager
//
// Class Description:
//
// Manager base class to handle event store and retrieve operation.
// Actual persistency implementation should be handled with derived classes.
//
// Each persistency package should implement derived classes of
// G4VHepMCIO, G4VMCTruthIO, G4VPHitIO, G4VPDigitIO, G4VPEventIO.
// Concreate G4PersistencyManager should implement the methods
// HepMCIO(), MCTruthIO(), HitIO(), DigitIO() and EventIO() to
// return the pointers of the above classes.
// G4PersistencyManager handles the sequence of the storing and
// retrieving of the persistent object of each type, along with
// the transaction handling.
//
// Retrieving a HepMC event:
//
//        G4PersistencyManager::Retrieve( HepMC::GenEvent*& )
//         |
//         |  ... StartRead() ...
//         |
//         |  ... Commit() ...
//         V
//
// Storing a Geant4 event:
//
//        G4PersistencyManager::Store( G4Pevent* )
//         |
//         |  ... StartUpdate() ...
//         |
//         |  ... MCTruthIO()->Store( MCTruth event ) ...
//         |
//         |  ... HitIO()->Store( hit_collection_of_event ) ...
//         |
//         |  ... DigitIO()->Store( digit_collection_of_event ) ...
//         |
//         |  ... EventIO()->Store( event with hits and digits ) ...
//         |
//         |  ... Commit() ...
//         V
//
// Retrieving a Geant4 event:
//
//        G4PersistencyManager::Retrieve( event )
//         |
//         |  ... StartRead() ...
//         |
//         |  ... EventIO()->Retrieve( event ) ...
//         |
//         |  ... Commit() ...
//         V
//
// Hit collection and digit collection of each detector component
// should be handled by detector specific I/O manager, which
// should be registered to the G4PersistencyCenter with
// AddHCIOmanager() and AddDCIOmanager().  Usually this is done
// through a command
//
//      /Persistency/Store/Using/HitIO <detector_io_manager_name>
//
// which is handled by G4PersistencyCenterMessenger.
//
// A static template declaration of G4HCIOentryT<class> must be
// implementated for each I/O manager.

// Author: Youhei Morita, 17.07.2001
// --------------------------------------------------------------------
#ifndef G4PERSISTENCYMANAGER_HH
#define G4PERSISTENCYMANAGER_HH 1

#include "G4Event.hh"

#include "G4VMCTruthIO.hh"
#include "G4HCIOcatalog.hh"
#include "G4DCIOcatalog.hh"
#include "G4VPEventIO.hh"
#include "G4VPHitIO.hh"
#include "G4VPDigitIO.hh"
#include "G4VTransactionManager.hh"
#include "G4VPersistencyManager.hh"

class G4PersistencyCenter;

class G4PersistencyManager : public G4VPersistencyManager
{
  friend class G4PersistencyCenter;

  public:

    G4PersistencyManager(G4PersistencyCenter* pc, const G4String& n);
      // Constructor

    virtual ~G4PersistencyManager();
      // Destructor

    virtual G4PersistencyManager* Create() { return nullptr; }
      // Create a new persistency manager.
      // To be used by G4PersistencyManagerT<>

    const G4String& GetName() { return nameMgr; }
      // Get the name of persistency manager

    virtual G4VPEventIO* EventIO() { return nullptr; }
      // Returns the current event I/O handling manager
      // Each derived class should return the pointer of actual manager

    virtual G4VPHitIO* HitIO() { return nullptr; }
      // Returns the current hit I/O handling manager
      // Each derived class should return the pointer of actual manager

    virtual G4VPDigitIO* DigitIO() { return nullptr; }
      // Returns the current digit I/O handling manager
      // Each derived class should return the pointer of actual manager

    virtual G4VMCTruthIO* MCTruthIO() { return nullptr; }
      // Returns the current MCTruth I/O handling manager
      // Each derived class should return the pointer of actual manager

    virtual G4VTransactionManager* TransactionManager() { return nullptr; }
      // Returns the current transaction manager
      // Each derived class should return the pointer of actual manager

    virtual void Initialize() {}
      // Initialize the persistency package.
      // Each derived class should implement the acutal initialization sequence

    void SetVerboseLevel(G4int v);
      // Set verbose level

    G4bool Store(const G4Event* evt);
      // Store the G4Event and its associated objects

    G4bool Retrieve(G4Event*& evt);
      // Retrieve the G4Event and its associated objects

    G4bool Store(const G4Run*) { return false; }
      // Not used

    G4bool Retrieve(G4Run*&) { return false; }
      // Not used

    G4bool Store(const G4VPhysicalVolume*) { return false; }
      // Not used

    G4bool Retrieve(G4VPhysicalVolume*&) { return false; }
      // Not used

  protected:

    static G4PersistencyManager* GetPersistencyManager();
      // Get the instance of persistency manager

  protected:

    G4PersistencyCenter* f_pc = nullptr;
    G4int m_verbose = 0;

  private:

    G4String nameMgr;
    G4bool f_is_initialized = false;
};

#endif
