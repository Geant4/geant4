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
// Rich advanced example for Geant4
// RichTbSD.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include "RichTbSD.hh"
#include "RichTbHit.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"

RichTbSD::RichTbSD(G4String DetName, G4int Numhpd , G4int NumSect,
                   G4int NumPixel,
                   G4String ) // CollectionName)
  :G4VSensitiveDetector(DetName),NumberOfSensitiveHpds(Numhpd),
   NumberOfSensitiveSectorsInHpd(NumSect),
   NumberOfSensitivePixelsInSector(NumPixel),
   HpdSDID(Numhpd*NumSect*NumPixel),HCID(-1)
{
  G4String HCname;
  collectionName.insert(HCname="RichTbHitsCollection");
     
}

RichTbSD::~RichTbSD(){

}

void RichTbSD::Initialize(G4HCofThisEvent*)
{
  RichTbHitCollection = new RichTbHitsCollection
                          (SensitiveDetectorName,collectionName[0]); 
  for (G4int jhpd =0; jhpd < NumberOfHpds; jhpd++){
    for(G4int jsect=0; jsect < NumberOfSiDetSectors; jsect++){
    for (G4int jpixel = 0; jpixel < NumberOfPadHpdSiPixels ; jpixel++ ){

        HpdSDID[NumberOfSensitiveSectorsInHpd*
                      NumberOfSensitivePixelsInSector*jhpd+
                      NumberOfSensitivePixelsInSector*jsect+jpixel]=-1;
    }
    }
  }

}


G4bool RichTbSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{
 


  G4Track* aTrack = aStep->GetTrack();
  if(aTrack ->GetDefinition()->GetPDGCharge() == 0.0 ) return false;
  G4double edep = aStep->GetTotalEnergyDeposit();

  if(edep < 0.0001) return false;
  G4StepPoint* PrePosition =  aStep->GetPreStepPoint();


  G4int CurrentHpdNumber =  PrePosition->GetTouchableHandle() 
    -> GetReplicaNumber(1);
  G4int CurrentSectorNumber= PrePosition->GetTouchableHandle()
    -> GetReplicaNumber(); 
  G4VPhysicalVolume* ROphysVol = ROhist -> GetVolume();
  G4int CurrentPixelNumber =  ROphysVol->GetCopyNo();
 
  G4ThreeVector HitAtPhotoCathode;
  if(aTrack->GetCreatorProcess() -> GetProcessName() == "PadHpdPhot") {
   HitAtPhotoCathode=aTrack->GetVertexPosition();
  }

  G4int CopyId= NumberOfSensitiveSectorsInHpd*
                      NumberOfSensitivePixelsInSector*CurrentHpdNumber+
                      NumberOfSensitivePixelsInSector*CurrentSectorNumber+
    CurrentPixelNumber;

  
  if( HpdSDID[CopyId ] == -1 ) {
  RichTbHit* newHit = new RichTbHit();
  newHit->SetEdep( edep );
  newHit->SetPos( aStep->GetPreStepPoint()->GetPosition() );
  newHit->SetPosPC(HitAtPhotoCathode);
  newHit->SetCurHpdNum ( CurrentHpdNumber );
  newHit->SetCurSectNum ( CurrentSectorNumber );
  newHit->SetCurPixNum (  CurrentPixelNumber );
  //For now the same SD detector cell is allowed to have
  //multiple hits from the same event. This is ok for
  // analogue readout for photoelectron counting.
  //This is to be changed in the future.
  G4int NumHits = RichTbHitCollection->insert( newHit );

  HpdSDID[CopyId]= NumHits -1 ;

  } else {
    // the current pixel is already hit in this event.
    // here we can add extra energy (adc counts) to the 
    // existing hit. But this is not relevant to the
    // current Tb analysis.
     (*RichTbHitCollection)[HpdSDID[CopyId]]->AddEdep( edep );
    
    G4cout << " Multiple hits in Hpd sector pixel " << CurrentHpdNumber  
           <<"   "<< CurrentSectorNumber<<"   "<< CurrentPixelNumber<< G4endl;
  }

  return true;
}

void  RichTbSD::EndOfEvent(G4HCofThisEvent*HCE){
  if( HCID < 0 ){
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
    HCE->AddHitsCollection( HCID, RichTbHitCollection );    

}

void  RichTbSD::clear(){} 

void   RichTbSD::DrawAll(){ } 

void   RichTbSD::PrintAll(){ } 




