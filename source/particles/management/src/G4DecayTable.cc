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
// $Id: G4DecayTable.cc 95470 2016-02-12 09:02:26Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
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
    G4cout << " G4DecayTable::Insert :: bad   G4VDecayChannel (mismatch parent) "
           << "       " << parent->GetParticleName()
           << " input:" << aChannel->GetParent()->GetParticleName() << G4endl;
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

G4VDecayChannel *G4DecayTable::SelectADecayChannel(G4double parentMass)
{
  // check if contents exist
  if (channels->size()<1) return 0;

  if(parentMass<0.) parentMass=parent->GetPDGMass(); 
  const size_t MAX_LOOP = 10000;
  for (size_t loop_counter=0; loop_counter <MAX_LOOP; ++loop_counter){
    G4double sumBR = 0.0;
    G4double r= G4UniformRand();
    // select decay channel
    G4VDecayChannelVector::iterator i;
    for (i = channels->begin(); i!= channels->end(); ++i) {
      sumBR += (*i)->GetBR();
      if ( !((*i)->IsOKWithParentMass(parentMass)) ) continue;
      if (r < sumBR) return (*i);
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











