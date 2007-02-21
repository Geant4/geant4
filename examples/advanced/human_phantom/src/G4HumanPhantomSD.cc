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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//
#include "G4HumanPhantomSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

G4HumanPhantomSD::G4HumanPhantomSD(G4String name)
  :G4VSensitiveDetector(name)
{
 G4String HCname;
 collectionName.insert(HCname="HumanPhantomCollection");
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

  //  G4cout <<bodypartName <<":" << edep/MeV<< G4endl; 

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
