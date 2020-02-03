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

G4DecayTable::G4DecayTable()
  :parent(nullptr), channels(nullptr)
{
  channels =  new G4VDecayChannelVector;
}

G4DecayTable::~G4DecayTable()
{
  // remove and delete all contents  
  G4VDecayChannelVector::iterator iCh;
  for (iCh = channels->begin(); iCh!= channels->end(); ++iCh) {
    delete (*iCh);
  }
  channels->clear();
  delete  channels;
  channels = nullptr;
  parent = nullptr;
}    

void G4DecayTable::Insert( G4VDecayChannel * aChannel){
  if (parent == nullptr) { 
    parent = (G4ParticleDefinition*)(aChannel->GetParent()); 
  }
  if (parent != aChannel->GetParent()) {
#ifdef G4VERBOSE
    G4cout << " G4DecayTable::Insert :: bad G4VDecayChannel (mismatch parent) "
           << "       " << parent->GetParticleName()
           << " input:" << aChannel->GetParent()->GetParticleName() << G4endl;
#endif
  } else {
    G4double br = aChannel->GetBR();
    G4VDecayChannelVector::iterator iCh;
    for (iCh = channels->begin(); iCh!= channels->end(); ++iCh) {
      if (br > (*iCh)->GetBR()) {
	channels->insert(iCh,aChannel);
	return;
      }
    }
    channels->push_back(aChannel);
  }
}

G4VDecayChannel *G4DecayTable::SelectADecayChannel(G4double parentMass)
{
  // check if contents exist
  if (channels->size()<1) return nullptr;

  if(parentMass < 0.) parentMass=parent->GetPDGMass(); 

  G4VDecayChannelVector::iterator iCh;
  G4double sumBR = 0.;
  for (iCh = channels->begin(); iCh!= channels->end(); ++iCh) {
    if ( !((*iCh)->IsOKWithParentMass(parentMass)) ) continue;
    sumBR += (*iCh)->GetBR();
  }
  if (sumBR <= 0.0) {
#ifdef G4VERBOSE
    G4cout << " G4DecayTable::SelectADecayChannel :: no possible DecayChannel"
           << "       " << parent->GetParticleName() << G4endl;
#endif
    return nullptr;
  }

  const size_t MAX_LOOP = 10000;
  for (size_t loop_counter=0; loop_counter <MAX_LOOP; ++loop_counter){
    G4double sum = 0.0;
    G4double br= sumBR*G4UniformRand();
    // select decay channel
    for (iCh = channels->begin(); iCh!= channels->end(); ++iCh) {
      sum += (*iCh)->GetBR();
      if ( !((*iCh)->IsOKWithParentMass(parentMass)) ) continue;
      if (br < sum) return (*iCh);
    }
  }
  return nullptr;
}

void G4DecayTable::DumpInfo() const
{
  G4cout << "G4DecayTable:  " << parent->GetParticleName() << G4endl;
  G4int index =0;
  G4VDecayChannelVector::iterator iCh;
  for (iCh = channels->begin(); iCh!= channels->end(); ++iCh) {
    G4cout << index << ": ";
    (*iCh)->DumpInfo();
    index += 1;
  }
  G4cout << G4endl;
}











