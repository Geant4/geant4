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
// $Id: G4SDManager.cc 96455 2016-04-15 08:44:06Z gcosmo $
//

#include "G4SDManager.hh"
#include "G4SDmessenger.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4VSensitiveDetector.hh"
#include "G4ios.hh"


G4ThreadLocal G4SDManager* G4SDManager::fSDManager = nullptr;

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
  delete HCtable;
  delete treeTop;
  theMessenger = nullptr;
  HCtable = nullptr;
  treeTop = nullptr;
  fSDManager = nullptr;
}

void G4SDManager::AddNewDetector(G4VSensitiveDetector*aSD)
{
  G4int numberOfCollections = aSD->GetNumberOfCollections();
  G4String pathName = aSD->GetPathName();
  if( pathName(0) != '/' ) pathName.prepend("/");
  if( pathName(pathName.length()-1) != '/' ) pathName += "/";
  treeTop->AddNewDetector(aSD,pathName);
  if(numberOfCollections<1) return;
  for(G4int i=0;i<numberOfCollections;i++)
  {
    G4String SDname = aSD->GetName();
    G4String DCname = aSD->GetCollectionName(i);
    AddNewCollection(SDname,DCname);
  }
  if( verboseLevel > 0 )
  {
    G4cout << "New sensitive detector <" << aSD->GetName()
         << "> is registered at " << pathName << G4endl;
  }
}

void G4SDManager::AddNewCollection(G4String SDname,G4String DCname)
{
  G4int i = HCtable->Registor(SDname,DCname);
  if(verboseLevel>0)
  {
    if(i<0) {
       if(verboseLevel>1) G4cout << "G4SDManager::AddNewCollection : the collection <"
        << SDname << "/" << DCname << "> has already been reginstered." << G4endl;
    }
    else
    {
      G4cout << "G4SDManager::AddNewCollection : the collection <"
       << SDname << "/" << DCname << "> is registered at " << i << G4endl;
    }
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

G4VSensitiveDetector* G4SDManager::FindSensitiveDetector(G4String dName, G4bool warning)
{
  G4String pathName = dName;
  if( pathName(0) != '/' ) pathName.prepend("/");
  return treeTop->FindSensitiveDetector(pathName, warning);
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


