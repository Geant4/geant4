// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DecayTable.cc,v 1.4 1999-05-25 02:01:08 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      27 July 1996 H.Kurashige
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
  channels->clearAndDestroy();
  delete channels;
}    

void G4DecayTable::Insert( G4VDecayChannel * aChannel){
  if (parent == 0) { parent = (G4ParticleDefinition*)(aChannel->GetParent()); }
  if (parent != aChannel->GetParent()) {
#ifdef G4VERBOSE
    G4cout << " G4DecayTable::Insert :: bad   G4VDecayChannel (mismatch parent) ";
    G4cout << "       " << parent->GetParticleName();
    G4cout << " input:" << aChannel->GetParent()->GetParticleName() << endl;
#endif
  } else {
    channels->insert(aChannel);
  }
}

G4VDecayChannel *G4DecayTable::SelectADecayChannel()
{
  // make cumlative table for branching ratio
  G4double sumBR = 0.0;
  G4int index;
  G4int    numberofchannels = channels->entries();
  G4double *cumBR = new G4double[numberofchannels];
  for (index=numberofchannels-1; index >=0 ; index--) {
    sumBR += ((*channels)(index))->GetBR();
    cumBR[index] = sumBR;
  }

  // select decay channel
  G4double r = sumBR* G4UniformRand();
  G4VDecayChannel *channel;
  for (index= numberofchannels-1; index >=0 ; index--) {
    if (r < cumBR[index]) {
      channel = (*channels)(index);
      break;
    }
  }
  delete [] cumBR;
  return channel;
}

void G4DecayTable::DumpInfo() const
{
  G4cout << "G4DecayTable:  " << parent->GetParticleName() << endl;
  G4int numberofchannels = channels->entries();
  for (G4int index= numberofchannels -1; index >=0 ; index -=1)
  {
    G4cout << index << ": ";
    ((*channels)(index))->DumpInfo();
  }
  G4cout << endl;
}











