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
// Med. Phys. 45  (2018) e722-e739
// Phys. Med. 31  (2015) 861-874
// Med. Phys. 37  (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157\u2013178
//
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

  // PROCESSING OF PROXIMITY FUNCTION t(x)
  
  // *************************************
  // Please select herebelow :
  // - the minimum value of x radius
  // - the maximum value of x radius
  // - the number of steps in radius
  
  G4double minRadius = 0.1 * nm;

  G4double maxRadius = 10000 * nm;
  G4int nRadiusSteps = 101;

  //
  
  auto analysisManager = G4AnalysisManager::Instance();

  G4double radius(minRadius);
  G4double stpRadius(std::pow(maxRadius/radius, 1./static_cast<G4double>(nRadiusSteps-1)));
  G4int step(nRadiusSteps);
  G4int noRadius(0);
   
  // 1) loop on radius
  
  while (step>0)
  {
    step--;
    noRadius=nRadiusSteps-step;
    
    //G4cout << "---radius/nm=" << radius/nm  << G4endl;

    // Computation of t(x)

    G4double tNum = 0.;
    G4double tDenom = 0.;
    G4int nbEdep = 0;

    // 2) loop on hits
    
    for ( G4int k=0; k<nofHits; k++ ) 
    { 
  
      G4ThreeVector hitPos =  (*fHitsCollection)[k]->GetPos();
      G4double hitNrj = (*fHitsCollection)[k]->GetEdep();

      //G4cout << "======> hit position x/nm =" << hitPos.x()/nm << G4endl; 
      //G4cout << "======> hit position y/nm =" << hitPos.y()/nm << G4endl; 
      //G4cout << "======> hit position z/nm =" << hitPos.z()/nm << G4endl; 
  
      // 3) loop on all other hits located within shell
        
      G4double localSum = 0.;

      for ( G4int i=0; i<nofHits; i++ ) 
      { 
        
        if ((*fHitsCollection)[i]->GetIncidentEnergy()>0)
          Einc = (*fHitsCollection)[i]->GetIncidentEnergy();

        G4ThreeVector localPosi = (*fHitsCollection)[i]->GetPos();
    
        if 
        (
         ( 
          (localPosi.x()-hitPos.x()) * (localPosi.x()-hitPos.x()) +
          (localPosi.y()-hitPos.y()) * (localPosi.y()-hitPos.y()) +
          (localPosi.z()-hitPos.z()) * (localPosi.z()-hitPos.z()) 
          < radius*stpRadius*radius*stpRadius
         )
         && 
         (
          (localPosi.x()-hitPos.x()) * (localPosi.x()-hitPos.x()) +
          (localPosi.y()-hitPos.y()) * (localPosi.y()-hitPos.y()) +
          (localPosi.z()-hitPos.z()) * (localPosi.z()-hitPos.z()) 
          >= radius*radius
         ) 
        )
           
        { 
          localSum = localSum + (*fHitsCollection)[i]->GetEdep() ;
          nbEdep = nbEdep + 1;
        }
      }

      tNum = tNum + localSum*hitNrj;
      tDenom = tDenom + hitNrj;
     

    } // loop on hits
  
    // fill ntuple including weighting
    // does not work with ntuple merging...

    analysisManager->FillNtupleDColumn(0,0, radius/nm);
    analysisManager->FillNtupleIColumn(0,1, noRadius);
    analysisManager->FillNtupleDColumn(0,2, nofHits);
    analysisManager->FillNtupleDColumn(0,3, nbEdep);
    analysisManager->FillNtupleDColumn(0,4, (tNum/tDenom)/eV);
    analysisManager->FillNtupleDColumn(0,5, (stpRadius*radius)/nm);
    analysisManager->FillNtupleDColumn(0,6, Einc/eV);
    analysisManager->AddNtupleRow();

    //G4cout << "---radius/nm=" << radius/nm << G4endl;
    //G4cout << "----end of radius--- " << radius/nm << G4endl;

    radius*=stpRadius;

  } // loop on radii       

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
