// $Id: G4RootEvent.cc,v 1.1 2002-12-04 02:44:29 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4RootEvent.cc
//
// History:
//   '02.5.7  Youhei Morita  Initial creation

#include "G4RootEvent.hh"

// Addtional Include:
#include "G4Event.hh"
#include "CLHEP/HepMC/GenEvent.h"
#include "G4MCTEvent.hh"
#include "G4Pevent.hh"

ClassImp(G4RootEvent)

// Implementation of Destructor #1
G4RootEvent::~G4RootEvent()
{}

// Implementation of MakeTransientObject
G4Pevent* G4RootEvent::MakeTransientObject()
{
  HepMC::GenEvent* hepevt = 0;
  G4MCTEvent*        mctevt = 0;
  G4Event*         g4evt  = 0;
  // G4Pevent* anEvt = new G4Pevent(hepevt, mctevt, g4evt);
  G4Pevent* anEvt = new G4Pevent(hepevt, g4evt);
  return anEvt;
}

// End of G4RootEvent.cc

