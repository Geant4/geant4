// $Id: G4PersistencyManager.hh,v 1.5 2002-12-04 13:57:29 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4PersistencyManager.hh
//
// History:
//   01.07.17  Youhei Morita  Initial creation (with "fadsclass")

#ifndef PERSISTENCY_MANAGER_HH
#define PERSISTENCY_MANAGER_HH 1

#include "G4Event.hh"
#ifndef WIN32
  #include "CLHEP/HepMC/GenEvent.h"
  #include "G4VHepMCIO.hh"
  #include "G4VMCTruthIO.hh"
#endif
#include "G4HCIOcatalog.hh"
#include "G4DCIOcatalog.hh"
#include "G4VPEventIO.hh"
#include "G4VPHitIO.hh"
#include "G4VPDigitIO.hh"
#include "G4VTransactionManager.hh"
#include <string>

class G4PersistencyCenter;

// Class inherited:
#include "G4VPersistencyManager.hh"

// Class Description:
//   Manager base class to handle event store and retrieve operation.
//   Actual persistency implementation should be handled with 
//   derived classes.
// 
//   Each persistency package should implement derived classes of
//   G4VHepMCIO, G4VMCTruthIO, G4VPHitIO, G4VPDigitIO, G4VPEventIO.
//   Concreate G4PersistencyManager should implement the methods
//   HepMCIO(), MCTruthIO(), HitIO(), DigitIO() and EventIO() to
//   return the pointers of the above classes.
//   G4PersistencyManager handles the sequence of the storing and
//   retrieving of the persistent object of each type, along with
//   the transaction handling.
// 
//   Retrieving a HepMC event:
// 
//        G4PersistencyManager::Retrieve( HepMC::GenEvent*& )
//         |
//         |  ... StartRead() ...
//         |
//         |  ... HepMCIO()->Retrieve() ...
//         |
//         |  ... Commit() ...
//         V
// 
//   Storing a Geant4 event:
// 
//        G4PersistencyManager::Store( G4Pevent* )
//         |
//         |  ... StartUpdate() ...
//         |
//         |  ... HepMCIO()->Store( HepMC event ) ...
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
//   Retrieving a Geant event:
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
//   Hit collection and digit collection of each detector component
//   should be handled by detector specific I/O manager, which
//   should be registered to the G4PersistencyCenter with
//   AddHCIOmanager() and AddDCIOmanager().  Usually this is done
//   through a command
// 
//      /Persistency/Store/Using/HitIO <detector_io_manager_name>
// 
//   which is handled by G4PersistencyCenterMessenger.
// 
//   A static template declaration of G4HCIOentryT<class> must be
//   implementated for each I/O manager.

class G4PersistencyManager
 : public G4VPersistencyManager
{
    friend class G4PersistencyCenter;

    public: // With description
      G4PersistencyManager(G4PersistencyCenter* pc, std::string n);
      // Constructor

      virtual ~G4PersistencyManager();
      // Destructor

    public: // With description
      virtual G4PersistencyManager* Create() {return 0;};
      // Create a new persistency manager.  To be used by G4PersistencyManagerT<>.

      std::string GetName() {return nameMgr;};
      // Get the name of persistency manager

      virtual G4VPEventIO* EventIO() { return 0; };
      // Returns the current event I/O handling manager
      // Each derived class should return the pointer of actual manager.

      virtual G4VPHitIO* HitIO() { return 0; };
      // Returns the current hit I/O handling manager
      // Each derived class should return the pointer of actual manager.

      virtual G4VPDigitIO* DigitIO() { return 0; };
      // Returns the current digit I/O handling manager
      // Each derived class should return the pointer of actual manager.
#ifndef WIN32
      virtual G4VHepMCIO* HepMCIO() { return 0; };
      // Returns the current HepMC I/O handling manager
      // Each derived class should return the pointer of actual manager.

      virtual G4VMCTruthIO* MCTruthIO() { return 0; };
      // Returns the current MCTruth I/O handling manager
      // Each derived class should return the pointer of actual manager.
#endif
      virtual G4VTransactionManager* TransactionManager() { return 0; };
      // Returns the current transaction manager
      // Each derived class should return the pointer of actual manager.

      virtual void Initialize() {};
      // Initialize the persistency package.
      // Each derived class should implement the acutal initialization sequence.

      void SetVerboseLevel(int v);
      // Set verbose level.

      G4bool Store(const G4Event* evt);
      // Store the G4Event and its associated objects

      G4bool Retrieve(G4Event*& evt);
      // Retrieve the G4Event and its associated objects
#ifndef WIN32
      G4bool Retrieve(HepMC::GenEvent*& evt, int id=-1);
      // retrieves HepMC GenEvent and its associated object.
      // To be used by generator/HepMCObjyReader.
#endif
      G4bool Store(const G4Run*) {return false;};
      // not used

      G4bool Retrieve(G4Run*&) {return false;};
      // not used

      G4bool Store(const G4VPhysicalVolume*) {return false;};
      // not used

      G4bool Retrieve(G4VPhysicalVolume*&) {return false;};
      // not used

    protected:
      static G4PersistencyManager* GetPersistencyManager();
      // Get the instance of persistency manager

    protected:
      G4PersistencyCenter* f_pc;
      int m_verbose;

    private:
      std::string      nameMgr;
      // GeneratorCenter* f_GenCenter;
      // G4MCTManager*      f_MCTman;
      G4bool             f_is_initialized;

}; // End of class G4PersistencyManager

#endif

