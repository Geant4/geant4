// $Id: G4PersistencyManager.cc,v 1.2 2002-11-29 09:10:55 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4PersistencyManager.cc
//
// History:
//   01.07.17  Youhei Morita  Initial creation (with "fadsclass")

#include "G4PersistencyManager.hh"

// Addtional Include:
#include <iomanip>
#include "G4PersistencyCenter.hh"

// Implementation of Constructor #1
G4PersistencyManager::G4PersistencyManager(G4PersistencyCenter* pc, std::string n)
 : f_pc(pc), nameMgr(n), f_is_initialized(false)
{
  m_verbose = f_pc->VerboseLevel();

  // std::cout << "G4PersistencyManager is constructed with name \"" << nameMgr
  //           << "\", " << this << ", verbose = " << m_verbose << std::endl;

  // f_GenCenter = GeneratorCenter::GetGeneratorCenter();
  // f_MCTman    = G4MCTManager::GetPointer();
}

// Implementation of Destructor #1
G4PersistencyManager::~G4PersistencyManager()
{}

// Implementation of GetG4PersistencyManager
G4PersistencyManager* G4PersistencyManager::GetG4PersistencyManager()
{
  return G4PersistencyCenter::GetG4PersistencyCenter()->CurrentG4PersistencyManager();
}

// Implementation of SetVerboseLevel
void G4PersistencyManager::SetVerboseLevel(int v)
{
  m_verbose = v;
  if ( m_verbose > 2 ) {
    std::cout << "G4PersistencyManager[\"" << nameMgr << "\"," << this
              << "]: verbose level is set to " << m_verbose << "."
              << std::endl;
  }
  if (   EventIO() != 0 )   EventIO()->SetVerboseLevel(m_verbose);
  if (   HepMCIO() != 0 )   HepMCIO()->SetVerboseLevel(m_verbose);
  if ( MCTruthIO() != 0 ) MCTruthIO()->SetVerboseLevel(m_verbose);
  if (     HitIO() != 0 )     HitIO()->SetVerboseLevel(m_verbose);
  if (   DigitIO() != 0 )   DigitIO()->SetVerboseLevel(m_verbose);
  if (TransactionManager() != 0) TransactionManager()->SetVerboseLevel(m_verbose);

  size_t i;

  G4HCIOcatalog* hcio = G4HCIOcatalog::GetG4HCIOcatalog();
  if ( hcio != 0 ) {
    hcio->SetVerboseLevel(m_verbose);
    for ( i = 0; i < hcio->NumberOfHCIOmanager(); i++ ) {
      hcio->GetHCIOmanager(i)->SetVerboseLevel(m_verbose);
    }
  }
  G4DCIOcatalog* dcio = G4DCIOcatalog::GetG4DCIOcatalog();
  if ( dcio != 0 ) {
    dcio->SetVerboseLevel(m_verbose);
    for ( i = 0; i < dcio->NumberOfDCIOmanager(); i++ ) {
      dcio->GetDCIOmanager(i)->SetVerboseLevel(m_verbose);
    }
  }
}

// Implementation of Store
G4bool G4PersistencyManager::Store(const G4Event* evt)
{
  if ( m_verbose > 2 ) {
    std::cout << "G4PersistencyManager::Store() is called for event# "
              << evt->GetEventID() << "." << std::endl;
  }

  bool is_store = f_pc->CurrentStoreMode("HepMC")   != kOff ||
                  f_pc->CurrentStoreMode("MCTruth") != kOff ||
                  f_pc->CurrentStoreMode("Hits")    != kOff ||
                  f_pc->CurrentStoreMode("Digits")  != kOff;

  if ( ! is_store ) return true;

  // Call package dependent Initialize()
  //
  if ( ! f_is_initialized ) {
    f_is_initialized = true;
    if ( m_verbose > 1 ) {
      std::cout << "G4PersistencyManager:: Initializing Transaction ... "
                << std::endl;
    }
    Initialize();
  }

  G4bool sthep = true, stmct = true;
  G4bool st1 = true, st2 = true, st3 = true;

  // Start event IO transaction
  //
  if ( TransactionManager()->StartUpdate() ) {
    if ( m_verbose > 2 ) {
      std::cout << "G4PersistencyManager: Update transaction started for event#"
                << evt->GetEventID() << "." << std::endl;
    }
  } else {
    std::cerr << "TransactionManager::Store(G4Event) - StartUpdate() failed."
              << std::endl;
    return false;
  }

  std::string file;

  // Store HepMC event
  //
  std::string obj = "HepMC";
  HepMC::GenEvent* hepevt = 0;
  if ( f_pc->CurrentStoreMode(obj) == kOn ) {

    //  Note: This part of code will not be activated until a method
    //  to obtain the current pointer of HepMC::GenEvent* become available.

    // if ( (hepevt = f_GenCenter->GetGenEvent()) !=0 ) {
    if ( hepevt !=0 ) {

      file = f_pc->CurrentWriteFile(obj);
      if ( TransactionManager()->SelectWriteFile(obj, file) ) {
        sthep = HepMCIO()->Store(hepevt);
        if ( sthep && m_verbose > 1 ) {
          std::cout << " -- File : " << file << " -- Event# "
                    << evt->GetEventID() << " -- HepMC Stored." << std::endl;
        }
      } else {
        sthep = false;
      }
    } // end of if ( hepevt != 0 )
  } else { // recycle or off
    // hepevt= f_GenCenter-> GetGenEvent();
  }

  // Store MCTruth event
  //
  obj = "MCTruth";
  G4MCTEvent* mctevt = 0;
  if ( f_pc->CurrentStoreMode(obj) == kOn ) {

    //  Note: This part of code will not be activated until a method
    //  to obtain the current pointer of G4MCTEvent* become available.

    // if ( (mctevt = f_MCTman->GetCurrentEvent()) != 0 ) {
    if ( mctevt != 0 ) {
      file = f_pc->CurrentWriteFile(obj);
      if ( TransactionManager()->SelectWriteFile(obj, file) ) {
        stmct = MCTruthIO()->Store(mctevt);
        if ( stmct && m_verbose > 1 ) {
          std::cout << " -- File : " << file << " -- Event# "
                  << evt->GetEventID() << " -- G4MCTEvent Stored." << std::endl;
        }
      } else {
        stmct = false;
      }
    } // end of if ( mctevt != 0 )
  }

  // Store hits collection
  //
  obj = "Hits";
  if ( f_pc->CurrentStoreMode(obj) == kOn ) {
    if ( G4HCofThisEvent* hc = evt->GetHCofThisEvent() ) {
      file = f_pc->CurrentWriteFile(obj);
      if ( TransactionManager()->SelectWriteFile(obj, file) ) {
        st1 = HitIO()->Store(hc);
        if ( st1 && m_verbose > 1 ) {
          std::cout << " -- File : " << file << " -- Event# "
                    << evt->GetEventID()
                    << " -- Hit Collections Stored." << std::endl;
        }
      } else {
        st1 = false;
      }
    }
  }

  // Store digits collection
  //
  obj = "Digits";
  if ( f_pc->CurrentStoreMode(obj) == kOn ) {
    if ( G4DCofThisEvent* dc = evt->GetDCofThisEvent() ) {
      file = f_pc->CurrentWriteFile(obj);
      if ( TransactionManager()->SelectWriteFile(obj, file) ) {
        st2 = DigitIO()->Store(dc);
        if ( st2 && m_verbose > 1 ) {
          std::cout << " -- File : " << file << " -- Event# "
                    << evt->GetEventID()
                    << " -- Digit Collections Stored." << std::endl;
        }
      } else {
        st2 = false;
      }
    }
  }

  // Store this G4EVENT
  //
  if ( hepevt!=0 || mctevt!=0 || evt!=0 ) {
    obj = "Hits";
    file = f_pc->CurrentWriteFile(obj);
    if ( TransactionManager()->SelectWriteFile(obj, file) ) {
      // st3 = EventIO()->Store(hepevt, mctevt, evt);
      st3 = EventIO()->Store(hepevt, evt);
      if ( st3 && m_verbose > 1 ) {
        std::cout << " -- File name: " << f_pc->CurrentWriteFile("Hits")
                  << " -- Event# "  << evt->GetEventID()
                  << " -- G4Pevent is Stored." << std::endl;
      }
    } else {
      st3 = false;
    }
  }

  G4bool st = sthep && stmct && st1 && st2 && st3;

  if ( st ) {
    TransactionManager()->Commit();
    if ( m_verbose > 0 )
      std::cout << "G4PersistencyManager: event# "
                << evt->GetEventID() << " is stored." << std::endl;
  } else {
    std::cerr << "G4PersistencyManager::Store(G4Event) - Transaction aborted."
              << std::endl;
    TransactionManager()->Abort();
  }

  return st;
}

// Implementation of Retrieve
G4bool G4PersistencyManager::Retrieve(G4Event*& evt)
{
  if ( m_verbose > 2 ) {
    std::cout << "G4PersistencyManager::Retrieve(G4Event*&) is called."
              << std::endl;
  }

  if ( f_pc->CurrentRetrieveMode("HepMC")   == false &&
       f_pc->CurrentRetrieveMode("MCTruth") == false &&
       f_pc->CurrentRetrieveMode("Hits")    == false &&
       f_pc->CurrentRetrieveMode("Digits")  == false ) {
    return true;
  }

  // Call package dependent Initialize()
  //
  if ( ! f_is_initialized ) {
    f_is_initialized = true;
    if ( m_verbose > 1 ) {
      std::cout << "G4PersistencyManager:: Initializing Transaction ... "
                << std::endl;
    }
    Initialize();
  }

  // Start event IO transaction
  //
  if ( TransactionManager()->StartRead() ) {
    if ( m_verbose > 2 ) {
      std::cout << "G4PersistencyManager: Read transaction started."
                << std::endl;
    }
  } else {
    std::cerr << "TransactionManager::Retrieve(G4Event) - StartRead() failed."
              << std::endl;
    return false;
  }

  G4bool st = false;
  std::string file;

  // Retrieve a G4EVENT
  //
  std::string obj = "Hits";
  if ( f_pc->CurrentRetrieveMode(obj) == true ) {
    file = f_pc->CurrentReadFile(obj);
    if ( TransactionManager()->SelectReadFile(obj, file) ) {
      st = EventIO()->Retrieve(evt);
      if ( st && m_verbose > 1 ) {
        std::cout << " -- File : " << file << " -- Event# "
                  << evt->GetEventID()
                  << " -- G4Event is Retrieved." << std::endl;
      }
    } else {
      st = false;
    }
  }

  if ( st ) {
    TransactionManager()->Commit();
  } else {
    std::cerr << "G4PersistencyManager::Retrieve() - Transaction aborted."
              << std::endl;
    TransactionManager()->Abort();
  }

  return st;
}

// Implementation of Retrieve
G4bool G4PersistencyManager::Retrieve(HepMC::GenEvent*& evt, int id)
{
  if ( m_verbose > 2 ) {
    std::cout << "G4PersistencyManager::Retrieve(HepMC::GenEvent*&) is called."
              << std::endl;
  }

  // Call package dependent Initialize()
  //
  if ( ! f_is_initialized ) {
    f_is_initialized = true;
    if ( m_verbose > 1 ) {
      std::cout << "G4PersistencyManager:: Initializing Transaction ... "
                << std::endl;
    }
    Initialize();
  }

  // Start event IO transaction
  //
  if ( TransactionManager()->StartRead() ) {
    if ( m_verbose > 2 ) {
      std::cout << "G4PersistencyManager: Read transaction started."
                << std::endl;
    }
  } else {
    std::cerr << "TransactionManager::Retrieve(HepMC) - StartRead() failed."
              << std::endl;
    return false;
  }

  G4bool st = false;
  std::string file;

  // Retrieve a HepMC GenEvent
  //
  std::string obj = "HepMC";
  if ( f_pc->CurrentRetrieveMode(obj) == true ) {
    file = f_pc->CurrentReadFile(obj);
    if ( TransactionManager()->SelectReadFile(obj, file) ) {
      st = HepMCIO()->Retrieve(evt, id);
      if ( st && m_verbose > 1 ) {
        std::cout << " -- File: " << file
                  << " - Event# " << HepMCIO()->LastEventID()
                  << " - HepMC event is Retrieved." << std::endl;
      }
    } else {
      st = false;
    }
  }

  if ( st ) {
    TransactionManager()->Commit();
  } else {
    std::cerr << "G4PersistencyManager::Retrieve(HepMC) - Transaction aborted."
              << std::endl;
    TransactionManager()->Abort();
  }

  return st;
}

// End of G4PersistencyManager.cc

