// $Id: G4MCTGenEvent.cc,v 1.2 2002-12-04 10:25:50 gcosmo Exp $
// ====================================================================
//
//   G4MCTGenEvent.cc
//
// ====================================================================

#ifndef WIN32

#include "G4MCTGenEvent.hh"

#include "CLHEP/HepMC/GenParticle.h"
#include "CLHEP/HepMC/GenVertex.h"

// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////////
G4MCTGenEvent::G4MCTGenEvent()
//////////////////////////
{
}

///////////////////////////
G4MCTGenEvent::~G4MCTGenEvent()
///////////////////////////
{
  eventList.clear();
}

/////////////////////////////////////////////////////////////
int G4MCTGenEvent::AddGenEvent(const HepMC::GenEvent* genevent)
/////////////////////////////////////////////////////////////
{
  eventList.push_back(const_cast<HepMC::GenEvent*>(genevent));
  return eventList.size();
}

/////////////////////////////////////
int G4MCTGenEvent::GetNofEvents() const
/////////////////////////////////////
{
  return eventList.size();
}

//////////////////////////////////////////////////////
const HepMC::GenEvent* G4MCTGenEvent::GetGenEvent(int i)
//////////////////////////////////////////////////////
{
  int size= eventList.size();
  if(i>=0 && i< size) return eventList[i];
  else return 0;
}


//////////////////////////////
void G4MCTGenEvent::ClearEvent()
//////////////////////////////
{
  eventList.clear();
}

////////////////////////////////////////////////////////////
void G4MCTGenEvent::Print(G4std::ostream& ostr) const
////////////////////////////////////////////////////////////
{
  int nev= eventList.size();
  for(int iev=0; iev<nev; iev++) {
    eventList[iev]-> print(ostr);
  }
}

#endif
