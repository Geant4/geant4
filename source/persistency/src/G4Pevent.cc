// $Id: G4Pevent.cc,v 1.2 2002-12-04 10:50:14 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4Pevent.cc
//
// History:
//   '01.11.18  Youhei Morita  Initial creation

#ifndef WIN32

#include "G4Pevent.hh"

//G4Pevent::G4Pevent( HepMC::GenEvent* hepevt, G4MCTEvent* mctevt, G4Event* g4evt )
// : f_hepevt(hepevt), genEventID(-1), f_mctevt(mctevt), f_g4evt(g4evt)

G4Pevent::G4Pevent( HepMC::GenEvent* hepevt, G4Event* g4evt )
 : f_hepevt(hepevt), genEventID(-1), f_g4evt(g4evt)
{
  m_id = g4evt->GetEventID();
  if (hepevt) genEventID= hepevt->event_number();
}

G4Pevent::~G4Pevent()
{
  delete f_g4evt;
}

#endif
