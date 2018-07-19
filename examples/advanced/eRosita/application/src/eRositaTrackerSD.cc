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
// $Id: eRositaTrackerSD.cc 107396 2017-11-10 08:28:08Z gcosmo $
//

#include "eRositaTrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4ParticleTypes.hh"

#include <iostream>
#include <fstream>

#include "AnalysisManager.hh"



eRositaTrackerSD::eRositaTrackerSD(G4String name)
  :G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="trackerCollection");
}


eRositaTrackerSD::~eRositaTrackerSD(){ }


void eRositaTrackerSD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new eRositaTrackerHitsCollection
    (SensitiveDetectorName,collectionName[0]); 
  static G4int HCID = -1;
  if (HCID < 0)
    { 
      HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); 
    }
  HCE->AddHitsCollection( HCID, trackerCollection ); 
}


G4bool eRositaTrackerSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  if(aStep->GetTrack()->GetDefinition() != G4Gamma::GammaDefinition()) return false;

//   G4double edep = aStep->GetTotalEnergyDeposit();
  G4double edep = aStep->GetPreStepPoint()->GetKineticEnergy();

  if (edep == 0.) return false;

  eRositaTrackerHit* newHit = new eRositaTrackerHit();
  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
  //newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchableHandle()
  //                                             ->GetCopyNumber());
  newHit->SetEdep(edep);
  newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
  trackerCollection->insert( newHit );
  
  //newHit->Print();
  //newHit->Draw();

  //ofstream out("ASCII");
  //newHit->PrintToFile(out);
  //out.close();

  return true;
}


void eRositaTrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  G4int NbHits = trackerCollection->entries();
 
  if (verboseLevel > 0) 
    { 
      
      G4cout << std::endl
	     << "Hits Collection: in this event they are " << NbHits 
	     << " hits in the tracker chambers: " << G4endl;
      for (G4int i=0;i<NbHits;i++) (*trackerCollection)[i]->Print();
      
    } 

  // ofstream out("ASCII");
  //if (!out.is_open())
  //  {
  //   G4cout <<"...opening ASCII file failed";
  // }
  
  double eTot = 0.;
  
  for (G4int i=0; i<NbHits; i++)
    {
      (*trackerCollection)[i]->PrintToFile();
      eTot += (*trackerCollection)[i]->GetEdep();
    };

//   if (eTot > 0.) AnalysisManager::Instance()->ScoreTot(eTot);

//out.close();
}




