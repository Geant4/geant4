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
//
// $Id: Tst50TrackerSD.cc,v 1.3 2003-01-16 14:11:51 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst50TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst50TrackerSD::Tst50TrackerSD(G4String name)
:G4VSensitiveDetector(name)
{ 
  G4String HCname;
  collectionName.insert(HCname="trackerCollection");
 hit_ID = new G4int[500];

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst50TrackerSD::~Tst50TrackerSD()
{
 delete [] hit_ID;
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst50TrackerSD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new Tst50TrackerHitsCollection
                          (SensitiveDetectorName,collectionName[0]); 
  for (G4int j=0;j<1;j++)
 {hit_ID [j]= -1;}; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Tst50TrackerSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
   if(edep==0.) return false;

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 
  theTouchable->MoveUpHistory();     
  
   
  if ( hit_ID[0]==-1)
   { 
     Tst50TrackerHit* newHit = new Tst50TrackerHit();
 
     newHit->SetEdep     (edep); 

     hit_ID[0] = trackerCollection->insert(newHit) - 1;
     // G4cout << " New Hit: "  << G4endl;
 }
  else
    { (*trackerCollection)[hit_ID[0]]->AddEnergy(edep); 
    //G4cout << " add Hit: "  << G4endl;
    }
 
 
 //newHit->Print();
  //newHit->Draw();

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst50TrackerSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection(HCID,trackerCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

