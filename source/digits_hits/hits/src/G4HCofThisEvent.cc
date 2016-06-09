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
// $Id: G4HCofThisEvent.cc,v 1.1 2003/10/03 10:18:00 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//

#include "G4HCofThisEvent.hh"

G4Allocator<G4HCofThisEvent> anHCoTHAllocator;

G4HCofThisEvent::G4HCofThisEvent()
{
  HC = new std::vector<G4VHitsCollection*>;
}

G4HCofThisEvent::G4HCofThisEvent(G4int cap)
{
  HC = new std::vector<G4VHitsCollection*>;
  for(G4int i=0;i<cap;i++)
  {
    HC->push_back((G4VHitsCollection*)0);
  }
}

G4HCofThisEvent::~G4HCofThisEvent()
{
  //HC->clearAndDestroy();
  for(size_t i=0;i<HC->size();i++)
  { delete (*HC)[i]; }
  HC->clear();
  delete HC;
}

void G4HCofThisEvent::AddHitsCollection(G4int HCID,G4VHitsCollection * aHC)
{
  if(HCID>=0 && HCID<G4int(HC->size()))
  { (*HC)[HCID] = aHC; }
}


