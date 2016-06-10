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
// $Id: G4DigiManager.cc 80987 2014-05-19 10:50:22Z gcosmo $
//

#include "G4DigiManager.hh"
#include "G4Event.hh"
#include "G4VHitsCollection.hh"
#include "G4VDigiCollection.hh"
#include "G4DMmessenger.hh"
#include "G4DCofThisEvent.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"


G4ThreadLocal G4DigiManager* G4DigiManager::fDManager = 0;

G4DigiManager* G4DigiManager::GetDMpointer()
{
  if(!fDManager)
  {
    fDManager = new G4DigiManager;
  }
  return fDManager;
}

G4DigiManager* G4DigiManager::GetDMpointerIfExist()
{ return fDManager; }

G4DigiManager::G4DigiManager():verboseLevel(0)
{ 
  theMessenger = new G4DMmessenger(this); 
  runManager = G4RunManager::GetRunManager();
  SDManager = G4SDManager::GetSDMpointer();
  DCtable = new G4DCtable;
}

G4DigiManager::~G4DigiManager()
{
  //DMtable.clearAndDestroy();
  for(G4int i=0;i<int(DMtable.size());i++)
  { delete DMtable[i]; }
  DMtable.clear();
  delete DCtable;
  delete theMessenger;
}

void G4DigiManager::AddNewModule(G4VDigitizerModule* DM)
{
  G4String DMname = DM->GetName();
  for(int j=0;j<int(DMtable.size());j++)
  {
    if(DMtable[j]==DM)
    { 
      G4cout << "<" << DMname << "> has already been registored." << G4endl; 
      return;
    }
  }
  if( verboseLevel > 0 )
  {
    G4cout << "New DigitizerModule <" << DMname
         << "> is registored." << G4endl;
  }
  DMtable.push_back(DM);

  G4int numberOfCollections = DM->GetNumberOfCollections();
  for(int i=0;i<numberOfCollections;i++)
  {
    G4String DCname = DM->GetCollectionName(i);
    if( DCtable->Registor(DMname,DCname) < 0 )
    { 
      G4cout << "DigiCollection <" << DCname 
           << "> has already been registored with "
           << DMname << " DigitizerModule." << G4endl;
    }
    else if( verboseLevel > 0 )
    {
      G4cout << "DigiCollection " << DCname 
           << " is registored. " << G4endl;
    }
  }
  
  runManager->SetDCtable(DCtable);
}

void G4DigiManager::Digitize(G4String mName)
{
  G4VDigitizerModule* aDM = FindDigitizerModule(mName);
  if(aDM)
  { aDM->Digitize(); }
  else
  { G4cout << "Unknown digitizer module <" << mName << ">. Digitize() ignored." << G4endl; }
}

G4VDigitizerModule* G4DigiManager::FindDigitizerModule(G4String mName)
{
  for(G4int i=0;i<int(DMtable.size());i++)
  {
    if(DMtable[i]->GetName() == mName) return DMtable[i];
  }
  return NULL;
}

const G4VHitsCollection* G4DigiManager::GetHitsCollection(G4int HCID,G4int eventID)
{
  const G4Event* evt = NULL;
  if(eventID==0)
  { evt = runManager->GetCurrentEvent(); }
  else
  { evt = runManager->GetPreviousEvent(eventID); }
  if(evt==NULL) return NULL;
  
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(HCE==NULL) return NULL;

  return HCE->GetHC(HCID);
}

const G4VDigiCollection* G4DigiManager::GetDigiCollection(G4int DCID,G4int eventID)
{
  const G4Event* evt = NULL;
  if(eventID==0)
  { evt = runManager->GetCurrentEvent(); }
  else
  { evt = runManager->GetPreviousEvent(eventID); }
  if(evt==NULL) return NULL;
  
  G4DCofThisEvent* DCE = evt->GetDCofThisEvent();
  if(DCE==NULL) return NULL;

  return DCE->GetDC(DCID);
}

G4int G4DigiManager::GetHitsCollectionID(G4String HCname)
{
  return SDManager->GetCollectionID(HCname);
}

G4int G4DigiManager::GetDigiCollectionID(G4String DCname)
{
  G4int i = DCtable->GetCollectionID(DCname);
  if(i==-2)
  { G4cout << "< " << DCname << "> is ambegious." << G4endl; }
  return i;
}

void G4DigiManager::SetDigiCollection(G4int DCID,G4VDigiCollection* aDC)
{
  const G4Event* consEvt = runManager->GetCurrentEvent();
  if(consEvt==NULL) 
  {
    G4cout << "G4DigiManager::SetDigiCollection --- "
         << "Event object is not available." << G4endl;
    return;
  }
  
  G4Event* evt = (G4Event*)consEvt;
  G4DCofThisEvent* DCE = evt->GetDCofThisEvent();
  if(DCE==NULL)
  {
    DCE = new G4DCofThisEvent(DCtable->entries());
    evt->SetDCofThisEvent(DCE);
    if(verboseLevel>0)
    { G4cout << "DCofThisEvent object is added to current G4Event." << G4endl; }
  }

  DCE->AddDigiCollection(DCID,aDC);

  if(verboseLevel>0)
  {
    G4cout << aDC->GetName() << " is stored at " << DCID 
         << "-th slot of G4DCofThisEvent." << G4endl;
  }
}
  
void G4DigiManager::SetVerboseLevel(G4int val)
{
  verboseLevel = val;
  for(G4int i=0;i<int(DMtable.size());i++)
  { DMtable[i]->SetVerboseLevel(val); }
}

void G4DigiManager::List() const
{
  for(G4int i=0;i<int(DMtable.size());i++)
  { G4cout << "   " << i << " : " << DMtable[i]->GetName() << G4endl; }
}



