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
//
//

#include "G4DCofThisEvent.hh"
#include <algorithm>

G4Allocator<G4DCofThisEvent>*& anDCoTHAllocator_G4MT_TLS_()
{
    G4ThreadLocalStatic G4Allocator<G4DCofThisEvent>* _instance = nullptr;
    return _instance;
}

G4DCofThisEvent::G4DCofThisEvent()
{ if (!anDCoTHAllocator_G4MT_TLS_()) anDCoTHAllocator_G4MT_TLS_() = new G4Allocator<G4DCofThisEvent>  ;
  DC = new std::vector<G4VDigiCollection*>;
}

G4DCofThisEvent::G4DCofThisEvent(G4int cap)
{ if (!anDCoTHAllocator_G4MT_TLS_()) anDCoTHAllocator_G4MT_TLS_() = new G4Allocator<G4DCofThisEvent>  ;
  DC = new std::vector<G4VDigiCollection*>;
  for(G4int i=0;i<cap;i++)
  {
    DC->push_back((G4VDigiCollection*)0);
  }
}

G4DCofThisEvent::~G4DCofThisEvent()
{ if (!anDCoTHAllocator_G4MT_TLS_()) anDCoTHAllocator_G4MT_TLS_() = new G4Allocator<G4DCofThisEvent>  ;
  //DC->clearAndDestroy();
  for(size_t i=0;i<DC->size();i++)
  { delete (*DC)[i]; }
  DC->clear();
  delete DC;
}

void G4DCofThisEvent::AddDigiCollection(G4int DCID,G4VDigiCollection * aDC)
{ if (!anDCoTHAllocator_G4MT_TLS_()) anDCoTHAllocator_G4MT_TLS_() = new G4Allocator<G4DCofThisEvent>  ;
  if(DCID>=0 && DCID<G4int(DC->size()))
  { (*DC)[DCID] = aDC; }
}

G4DCofThisEvent::G4DCofThisEvent(const G4DCofThisEvent& rhs)
{
    if ( !anDCoTHAllocator_G4MT_TLS_() ) anDCoTHAllocator_G4MT_TLS_() = new G4Allocator<G4DCofThisEvent>;
    DC = new std::vector<G4VDigiCollection*>(rhs.DC->size());
    for ( unsigned int i = 0 ; i < rhs.DC->size() ; ++i )
        *(DC->at(i)) = *(rhs.DC->at(i));
}

G4DCofThisEvent& G4DCofThisEvent::operator=(const G4DCofThisEvent& rhs)
{
    if ( this == &rhs ) return *this;
    if ( !anDCoTHAllocator_G4MT_TLS_() ) anDCoTHAllocator_G4MT_TLS_() = new G4Allocator<G4DCofThisEvent>;
    for ( std::vector<G4VDigiCollection*>::const_iterator it = DC->begin() ;
         it != DC->end() ; ++it )
    {
        delete *it;
    }
    DC->resize(rhs.DC->size());
    for ( unsigned int i = 0 ; i < rhs.DC->size() ; ++i )
        *(DC->at(i)) = *(rhs.DC->at(i));
    return *this;
}

