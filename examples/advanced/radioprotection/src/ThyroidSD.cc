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
//    ********************************
//    *                              *  
//    *         ThyroidSD.cc         *
//    *                              *
//    ********************************

#include "ThyroidSD.hh"
#include "ThyroidHit.hh"
#include "ThyroidDetectorConstruction.hh"

#include "G4Track.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"

//....

ThyroidSD::ThyroidSD(G4String name)
  :G4VSensitiveDetector(name)
{
 G4String HCname;
 collectionName.insert(HCname="RightThyroidHitsCollection");
 m_pRightThyroidHitsCollection = NULL;
}

//....
ThyroidSD:: ~ThyroidSD(){}

 void ThyroidSD::Initialize(G4HCofThisEvent*HCE)
{
 m_pRightThyroidHitsCollection = new
ThyroidHitsCollection(SensitiveDetectorName,collectionName[0]);  
 
 }

//....

G4bool ThyroidSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
 printf ("\n sono in SD \n");

 G4double edep = aStep->GetTotalEnergyDeposit();
 if(edep==0.)
	return false;
 G4VPhysicalVolume* physVol = ROhist->GetVolume();
 G4VPhysicalVolume* mothVol = ROhist->GetVolume(1);

	ThyroidHit* RightThyroidHit = new
ThyroidHit(physVol->GetLogicalVolume());       
 
   
	RightThyroidHit->SetEdep(edep);
        RightThyroidHit->SetPos(physVol->GetTranslation());
   


	return true;
 }

void ThyroidSD::EndOfEvent(G4HCofThisEvent*HCE)
{ static G4int HCID = -1;
  if(HCID<0)
  	{ 
	HCID = GetCollectionID(0); 
	}
  HCE->AddHitsCollection(HCID,m_pRightThyroidHitsCollection);
}

 

//....

void ThyroidSD::clear()
{
} 

//....

void ThyroidSD::DrawAll()
{
} 

//....

void ThyroidSD::PrintAll()
{
}







