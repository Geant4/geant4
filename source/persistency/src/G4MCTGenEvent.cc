//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
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
