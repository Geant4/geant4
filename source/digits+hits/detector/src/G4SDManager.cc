// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SDManager.cc,v 1.1.10.1.2.1 1999/12/08 17:33:16 gunter Exp $
// GEANT4 tag $Name: geant4-03-00 $
//

#include "G4SDManager.hh"
#include "G4SDmessenger.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4VSensitiveDetector.hh"
#include "G4ios.hh"


G4SDManager* G4SDManager::fSDManager = 0;

G4SDManager* G4SDManager::GetSDMpointer()
{
  if(!fSDManager)
  {
    fSDManager = new G4SDManager;
  }
  return fSDManager;
}

G4SDManager* G4SDManager::GetSDMpointerIfExist()
{ return fSDManager; }

G4SDManager::G4SDManager():verboseLevel(0)
{
  G4String topName = "/";
  treeTop = new G4SDStructure(topName);
  theMessenger = new G4SDmessenger(this);
  HCtable = new G4HCtable;
}

G4SDManager::~G4SDManager()
{
  delete theMessenger;
}

void G4SDManager::AddNewDetector(G4VSensitiveDetector*aSD)
{
  G4String pathName = aSD->GetPathName();
  G4int numberOfCollections = aSD->GetNumberOfCollections();
  if( pathName(0) != '/' ) pathName.prepend("/");
  if( pathName(pathName.length()-1) != '/' ) pathName += "/";
  treeTop->AddNewDetector(aSD,pathName);
  for(int i=0;i<numberOfCollections;i++)
  {
    G4String SDname = aSD->GetName();
    G4String DCname = aSD->GetCollectionName(i);
    HCtable->Registor(SDname,DCname);
  }
  if( verboseLevel > 0 )
  {
    G4cout << "New sensitive detector <" << aSD->GetName()
         << "> is registored at " << pathName << G4endl;
  }
}

G4HCofThisEvent* G4SDManager::PrepareNewEvent()
{
  G4HCofThisEvent* HCE = new G4HCofThisEvent(HCtable->entries());
  treeTop->Initialize(HCE);
  return HCE;
}

void G4SDManager::TerminateCurrentEvent(G4HCofThisEvent* HCE)
{
  treeTop->Terminate(HCE);
}

void G4SDManager::Activate(G4String dName, G4bool activeFlag)
{
  G4String pathName = dName;
  if( pathName(0) != '/' ) pathName.prepend("/");
  treeTop->Activate(pathName,activeFlag);
}

G4VSensitiveDetector* G4SDManager::FindSensitiveDetector(G4String dName)
{
  G4String pathName = dName;
  if( pathName(0) != '/' ) pathName.prepend("/");
  return treeTop->FindSensitiveDetector(pathName);
}

G4int G4SDManager::GetCollectionID(G4String colName)
{
  G4int id = HCtable->GetCollectionID(colName);
  if(id==-1)
  { G4cout << "<" << colName << "> is not found." << G4endl; }
  else if(id==-2)
  { G4cout << "<" << colName << "> is ambiguous." << G4endl; }
  return id;
}

G4int G4SDManager::GetCollectionID(G4VHitsCollection* aHC)
{
  G4String HCname = aHC->GetSDname();
  HCname += "/";
  HCname += aHC->GetName();
  return GetCollectionID(HCname);
}


