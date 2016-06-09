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
//
//
// $Id: G4DCofThisEvent.cc,v 1.3 2004/06/11 14:10:31 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//

#include "G4DCofThisEvent.hh"

G4Allocator<G4DCofThisEvent> anDCoTHAllocator;

G4DCofThisEvent::G4DCofThisEvent()
{
  DC = new std::vector<G4VDigiCollection*>;
}

G4DCofThisEvent::G4DCofThisEvent(G4int cap)
{
  DC = new std::vector<G4VDigiCollection*>;
  for(G4int i=0;i<cap;i++)
  {
    DC->push_back((G4VDigiCollection*)0);
  }
}

G4DCofThisEvent::~G4DCofThisEvent()
{
  //DC->clearAndDestroy();
  for(size_t i=0;i<DC->size();i++)
  { delete (*DC)[i]; }
  DC->clear();
  delete DC;
}

void G4DCofThisEvent::AddDigiCollection(G4int DCID,G4VDigiCollection * aDC)
{
  if(DCID>=0 && DCID<G4int(DC->size()))
  { (*DC)[DCID] = aDC; }
}

