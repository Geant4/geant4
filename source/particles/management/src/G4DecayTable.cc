// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DecayTable.cc,v 1.6 2000-02-25 07:36:23 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      27 July 1996 H.Kurashige
// ----------------------------------------
//      implementation for STL          14 Feb. 2000 H.Kurashige
// ------------------------------------------------------------

#include "globals.hh"
#include "G4DecayTable.hh"
#include "Randomize.hh"

G4DecayTable::G4DecayTable():parent(0)
{
  channels =  new G4VDecayChannelVector;
}

G4DecayTable::~G4DecayTable()
{
  // remove and delete all contents  
  G4VDecayChannelVector::iterator i;
  for (i = channels->begin(); i!= channels->end(); ++i) {
    delete (*i);
  }
  channels->clear();
  delete  channels;
  channels = 0;
}    

void G4DecayTable::Insert( G4VDecayChannel * aChannel){
  if (parent == 0) { parent = (G4ParticleDefinition*)(aChannel->GetParent()); }
  if (parent != aChannel->GetParent()) {
#ifdef G4VERBOSE
    G4cout << " G4DecayTable::Insert :: bad   G4VDecayChannel (mismatch parent) ";
    G4cout << "       " << parent->GetParticleName();
    G4cout << " input:" << aChannel->GetParent()->GetParticleName() << G4endl;
#endif
  } else {
    G4double r = aChannel->GetBR();
    G4VDecayChannelVector::iterator i;
    for (i = channels->begin(); i!= channels->end(); ++i) {
      if (r > (*i)->GetBR()) {
	channels->insert(i,aChannel);
	return;
      }
    }
    channels->push_back(aChannel);
  }
}

G4VDecayChannel *G4DecayTable::SelectADecayChannel()
{
  // check if contents exist
  if (channels->size()<1) return 0;

  while (1) {
    G4double sumBR = 0.0;
    G4double r= G4UniformRand();
    // select decay channel
    G4VDecayChannelVector::iterator i;
    for (i = channels->begin(); i!= channels->end(); ++i) {
      sumBR += (*i)->GetBR();
      if (r < sumBR) {
	return (*i);
      }
    }
  }
  return 0;
}

void G4DecayTable::DumpInfo() const
{
  G4cout << "G4DecayTable:  " << parent->GetParticleName() << G4endl;
  G4int index =0;
  G4VDecayChannelVector::iterator i;
  for (i = channels->begin(); i!= channels->end(); ++i) {
    G4cout << index << ": ";
    (*i)->DumpInfo();
    index +=1;
  }
  G4cout << G4endl;
}











