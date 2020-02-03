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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
//
#include "G4HumanPhantomSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

G4HumanPhantomSD::G4HumanPhantomSD(const G4String& name,
                         const G4String& hitsCollectionName)
  :G4VSensitiveDetector(name)
{
 G4String HCname;
 collectionName.insert(hitsCollectionName);
}

G4HumanPhantomSD::~G4HumanPhantomSD()
{ 
}

void G4HumanPhantomSD::Initialize(G4HCofThisEvent* HCE)
{
  collection = new G4HumanPhantomHitsCollection
                          (SensitiveDetectorName,collectionName[0]); 
  static G4int HCID = -1;
  if(HCID<0)
  { 
  HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, collection ); 
}

G4bool G4HumanPhantomSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{  
  G4double edep = aStep->GetTotalEnergyDeposit();

  if(edep==0.) return false;
 
  G4String bodypartName = aStep->GetPreStepPoint()->GetTouchable()
   			->GetVolume()->GetLogicalVolume()->GetName();

  // G4cout <<bodypartName <<":" << edep/MeV<< G4endl; 

  G4HumanPhantomHit* newHit = new G4HumanPhantomHit();
  newHit->SetEdep(edep);
  newHit->SetBodyPartID(bodypartName);
  collection->insert(newHit);
  
 return true;
}

void G4HumanPhantomSD::EndOfEvent(G4HCofThisEvent*)
{

// G4int NbHits = collection->entries();
//      G4cout << "\n-------->Hits Collection: in this event they are " << NbHits 
//             << " hits in the tracker chambers: " << G4endl;
//      for (G4int i=0;i<NbHits;i++) (*collection)[i]->Print();
}
