// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PersistentEventMan.cc,v 1.15 1999/12/02 19:40:51 morita Exp $
// GEANT4 tag $Name: geant4-03-01 $
//
// class G4PersistentEventMan 
//
// Implementation for concrete G4PersistentEventMan.
//
// History:
// 98.01.08 Y.Morita  Initial version
// 98.10.30 Y.Morita  Splitted from G4PersistencyManager

#include "G4PersistentEventMan.hh"

#include "G4PersistentHitMan.hh"
#include "G4PersistentDigitMan.hh"

#include "HepODBMS/clustering/HepDbApplication.h"

#include "G4PersistentHitMan.hh"
#include "G4PersistentDigitMan.hh"
#include "G4Event.hh"
#include "G4ios.hh"

G4PersistentEventMan::G4PersistentEventMan()
 : f_PHCMan(NULL), f_PDCMan(NULL), f_currentEventID(0)
{;}

G4PersistentEventMan::G4PersistentEventMan(
                            G4PersistentHitMan*   aPHCMan,
                            G4PersistentDigitMan* aPDCMan )
 : f_PHCMan(aPHCMan), f_PDCMan(aPDCMan)
{;}

G4PersistentEventMan::~G4PersistentEventMan()
{;}

//----------------------------------------------------------------------------

G4bool G4PersistentEventMan::Store( HepDbApplication* dbApp,
                                    const G4Event* anEvent )
{
  HepRef(G4PHCofThisEvent) aPHC = f_PHCMan->GetCurrentPHCofThisEvent();
  HepRef(G4PDCofThisEvent) aPDC = f_PDCMan->GetCurrentPDCofThisEvent();

  // Create persistent event
  f_currentPEvent = new(f_container) G4PEvent(anEvent, aPHC, aPDC);

  if( aPHC != NULL )
  { f_PHCMan->SetCurrentPHCofThisEvent(NULL); }

  if( aPDC != NULL )
  { f_PDCMan->SetCurrentPDCofThisEvent(NULL); }

  if( f_currentPEvent != NULL )
  {
    f_currentEventID = anEvent->GetEventID();
    return true;
  }
  else
  {
    f_currentEventID = -1;
    return false;
  }
}

G4bool G4PersistentEventMan::Retrieve( HepDbApplication* dbApp,
                                       G4Event*& anEvent )
{
  G4bool theStatus = false;
  anEvent = NULL;

  ooItr(G4PEvent) pevt_iterator;

  // set the new scan scope if the DB/Container has been changed by
  // G4PersistencyManager
  if( f_container != f_currentContainer )
  {
    f_currentContainer = f_container;
    pevt_iterator.scan(f_currentContainer);
  }

  // Retrieve "next" G4PEvent in this scope and make a G4Event
  if( pevt_iterator.next() )
  {
    G4Event* anEvt = pevt_iterator->MakeTransientObject();
    if( anEvt != NULL )
    {
      anEvent = anEvt;
      theStatus = true;
      f_currentEventID = anEvent->GetEventID();
    }
  }

  return theStatus;
}

