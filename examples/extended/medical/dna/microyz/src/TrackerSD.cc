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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publications:
// Phys. Med. 31 (2015) 861-874
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file TrackerSD.cc
/// \brief Implementation of the TrackerSD class

#include "Analysis.hh"
#include "TrackerSD.hh"
#include "Randomize.hh"
#include "G4SDManager.hh"

#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::TrackerSD(const G4String& name,
                     const G4String& hitsCollectionName) 
:G4VSensitiveDetector(name),
fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::~TrackerSD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection 
    = new TrackerHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  
  hce->AddHitsCollection( hcID, fHitsCollection ); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackerSD::ProcessHits(G4Step* aStep, 
                                     G4TouchableHistory*)
{  
  // energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();

  if (edep==0.) return false;

  TrackerHit* newHit = new TrackerHit();

  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
  newHit->SetEdep(edep);
  newHit->SetPos (aStep->GetPostStepPoint()->GetPosition());

  if (aStep->GetTrack()->GetTrackID()==1&&aStep->GetTrack()->GetParentID()==0)
    newHit->SetIncidentEnergy(aStep->GetTrack()->GetVertexKineticEnergy());

  fHitsCollection->insert( newHit );

  //newHit->Print();

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  G4int nofHits = fHitsCollection->entries();

  G4double Einc=0;
  
  /*
  G4cout << G4endl
  << "-------->Hits Collection: in this event they are " 
  << nofHits 
  << " hits in the target volume " << G4endl;
  */
  
  // PROCESSING OF MICRODOSIMETRY Y & Z SPECTRA
  
  // *************************************
  // Please select herebelow :
  // the radius of the target sphere:
  // variable name = radius
  // it is set to 5 nm by default)
  //
  
  G4double radius = 5*nm;

  //
 
  //***************
  // y and z
  //***************
  
  // select random hit
  G4int randHit=0; // Runs from 0 to number of hits - 1
  randHit = static_cast<G4int>( G4UniformRand()*nofHits );
  
  /*
  G4cout 
  << "======> random selection of hit number randHit =" 
  << randHit << G4endl;
  */
  
  // get selected random hit position
  G4ThreeVector hitPos =  (*fHitsCollection)[randHit]->GetPos();
//G4cout << "======> random hit position x/nm =" << hitPos.x()/nm << G4endl; 
//G4cout << "======> random hit position y/nm =" << hitPos.y()/nm << G4endl; 
//G4cout << "======> random hit position z/nm =" << hitPos.z()/nm << G4endl; 
  
  // set random position of center of sphere within radius
  G4double chord = 4.*radius/3;
  G4double density = 1 * g/cm3;
  G4double mass = (4./3)*CLHEP::pi*radius*radius*radius*density;
  
  // random placement of sphere: method 1
  /*  
  G4ThreeVector randDir = G4RandomDirection();
  G4double randRadius = G4UniformRand()*radius;
  G4ThreeVector randCenterPos = randRadius*randDir + hitPos;
  */  

  // random placement of sphere: method 2

  G4double xRand = 1.01*radius;
  G4double yRand = 1.01*radius;
  G4double zRand = 1.01*radius;
  G4double randRad = 1.01*radius;
  do
  {
    xRand = (2*G4UniformRand()-1)*radius;
    yRand = (2*G4UniformRand()-1)*radius;
    zRand = (2*G4UniformRand()-1)*radius;
    randRad = std::sqrt( xRand*xRand+yRand*yRand+zRand*zRand );
  }
  while (randRad>radius);

  G4ThreeVector 
    randCenterPos(xRand+hitPos.x(),yRand+hitPos.y(),zRand+hitPos.z());

  // search for neighbouring hits in the sphere and cumulate deposited energy 
  //  in epsilon
  G4double epsilon = 0;
  G4int nbEdep = 0;
  
  for ( G4int i=0; i<nofHits; i++ ) 
  { 

    if ((*fHitsCollection)[i]->GetIncidentEnergy()>0)
      Einc = (*fHitsCollection)[i]->GetIncidentEnergy();
    
    G4ThreeVector localPos = (*fHitsCollection)[i]->GetPos();
    
    // G4cout << i << " " << (*fHitsCollection)[i] << G4endl;
    // G4cout << i << " " << (*fHitsCollection)[i]->GetEdep()/eV << G4endl;
    
    if ( 
        (localPos.x()-randCenterPos.x()) * (localPos.x()-randCenterPos.x()) +
        (localPos.y()-randCenterPos.y()) * (localPos.y()-randCenterPos.y()) +
        (localPos.z()-randCenterPos.z()) * (localPos.z()-randCenterPos.z()) 
         <= radius*radius
       ) 
       
    { 
      epsilon = epsilon + (*fHitsCollection)[i]->GetEdep() ;
      nbEdep = nbEdep+1;
    }
       
  }

  // for testing only
  /*
  G4cout << "======> for hit number #" << randHit <<
  ", we collect " 
  << nbEdep << " energy depositions in a sphere of radius " 
  << radius/nm << " nm and mass " 
  << mass/kg << " kg for a total of " 
  << epsilon/eV << " eV or " 
  << (epsilon/joule)/(mass/kg) << " Gy" << G4endl;
  G4cout << "-" << G4endl;
  */
  
  /*
  FILE* myFile;
  myFile=fopen("yz.txt","a");
  fprintf(myFile,"%e %e %e\n",radius/nm,(epsilon/eV)/(chord/nm),
   (epsilon/joule)/(mass/kg));
  fclose(myFile);
  */

  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // fill ntuple including weighting
  analysisManager->FillNtupleDColumn(0, radius/nm);
  analysisManager->FillNtupleDColumn(2, nofHits);
  analysisManager->FillNtupleDColumn(3, nbEdep);
  analysisManager->FillNtupleDColumn(4, (epsilon/eV)/(chord/nm));
  analysisManager->FillNtupleDColumn(5, (epsilon/mass)/gray);
  analysisManager->FillNtupleDColumn(6, Einc/eV);
  analysisManager->AddNtupleRow();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
