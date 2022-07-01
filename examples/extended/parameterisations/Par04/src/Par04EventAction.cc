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
#include "Par04EventAction.hh"
#include <CLHEP/Units/SystemOfUnits.h>   // for GeV
#include <CLHEP/Vector/ThreeVector.h>    // for Hep3Vector
#include <stddef.h>                      // for size_t
#include <G4Exception.hh>                // for G4Exception, G4ExceptionDesc...
#include <G4ExceptionSeverity.hh>        // for FatalException
#include <G4GenericAnalysisManager.hh>   // for G4GenericAnalysisManager
#include <G4PrimaryParticle.hh>          // for G4PrimaryParticle
#include <G4PrimaryVertex.hh>            // for G4PrimaryVertex
#include <G4SystemOfUnits.hh>            // for GeV
#include <G4THitsCollection.hh>          // for G4THitsCollection
#include <G4ThreeVector.hh>              // for G4ThreeVector
#include <G4Timer.hh>                    // for G4Timer
#include <G4UserEventAction.hh>          // for G4UserEventAction
#include <algorithm>                     // for max
#include <ostream>                       // for basic_ostream::operator<<
#include "G4AnalysisManager.hh"          // for G4AnalysisManager
#include "G4Event.hh"                    // for G4Event
#include "G4EventManager.hh"             // for G4EventManager
#include "G4HCofThisEvent.hh"            // for G4HCofThisEvent
#include "G4SDManager.hh"                // for G4SDManager
#include "Par04DetectorConstruction.hh"  // for Par04DetectorConstruction
#include "Par04Hit.hh"                   // for Par04Hit, Par04HitsCollection
#include <cmath>                         // for log10

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04EventAction::Par04EventAction(Par04DetectorConstruction* aDetector)
  : G4UserEventAction()
  , fHitCollectionID(-1)
  , fTimer()
  , fDetector(aDetector)
{
  fCellNbRho = aDetector->GetMeshNbOfCells().x();
  fCellNbPhi = aDetector->GetMeshNbOfCells().y();
  fCellNbZ   = aDetector->GetMeshNbOfCells().z();
  fCalEdep.reserve(fCellNbRho * fCellNbPhi * fCellNbZ);
  fCalRho.reserve(fCellNbRho * fCellNbPhi * fCellNbZ);
  fCalPhi.reserve(fCellNbRho * fCellNbPhi * fCellNbZ);
  fCalZ.reserve(fCellNbRho * fCellNbPhi * fCellNbZ);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04EventAction::~Par04EventAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04EventAction::BeginOfEventAction(const G4Event*) { fTimer.Start(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04EventAction::EndOfEventAction(const G4Event* aEvent)
{
  fTimer.Stop();

  // Get hits collection ID (only once)
  if(fHitCollectionID == -1)
  {
    fHitCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("hits");
  }
  // Get hits collection
  auto hitsCollection =
    static_cast<Par04HitsCollection*>(aEvent->GetHCofThisEvent()->GetHC(fHitCollectionID));

  if(hitsCollection == nullptr)
  {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << fHitCollectionID;
    G4Exception("Par04EventAction::GetHitsCollection()", "MyCode0001", FatalException, msg);
  }

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  // Retrieve only once detector dimensions
  if(fCellSizeZ == 0)
  {
    fCellSizeZ   = fDetector->GetMeshSizeOfCells().z();
    fCellSizePhi = fDetector->GetMeshSizeOfCells().y();
    fCellSizeRho = fDetector->GetMeshSizeOfCells().x();
    fCellNbRho   = fDetector->GetMeshNbOfCells().x();
    fCellNbPhi   = fDetector->GetMeshNbOfCells().y();
    fCellNbZ     = fDetector->GetMeshNbOfCells().z();
  }

  // Retrieve information from primary vertex and primary particle
  // To calculate shower axis and entry point to the detector
  auto primaryVertex =
    G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetPrimaryVertex();
  auto primaryParticle   = primaryVertex->GetPrimary(0);
  G4double primaryEnergy = primaryParticle->GetTotalEnergy();
  // Estimate from vertex and particle direction the entry point to the detector
  // Calculate entrance point to the detector located at z = 0
  auto primaryDirection = primaryParticle->GetMomentumDirection();
  auto primaryEntrance =
    primaryVertex->GetPosition() - primaryVertex->GetPosition().z() * primaryDirection;

  // Resize back to initial mesh size
  fCalEdep.resize(fCellNbRho * fCellNbPhi * fCellNbZ);
  fCalRho.resize(fCellNbRho * fCellNbPhi * fCellNbZ);
  fCalPhi.resize(fCellNbRho * fCellNbPhi * fCellNbZ);
  fCalZ.resize(fCellNbRho * fCellNbPhi * fCellNbZ);

  // Fill histograms
  Par04Hit* hit                  = nullptr;
  G4double hitEn                 = 0;
  G4double totalEnergy           = 0;
  G4int hitZ                     = -1;
  G4int hitRho                   = -1;
  G4int hitPhi                   = -1;
  G4int hitType                  = -1;
  G4int numNonZeroThresholdCells = 0;
  G4double tDistance = 0., rDistance = 0., phiDistance = 0.;
  G4double tFirstMoment = 0., tSecondMoment = 0.;
  G4double rFirstMoment = 0., rSecondMoment = 0.;
  G4double phiMean = 0.;
  for(size_t iHit = 0; iHit < hitsCollection->entries(); iHit++)
  {
    hit     = static_cast<Par04Hit*>(hitsCollection->GetHit(iHit));
    hitZ    = hit->GetZid();
    hitRho  = hit->GetRhoId();
    hitPhi  = hit->GetPhiId();
    hitEn   = hit->GetEdep();
    hitType = hit->GetType();
    if(hitEn > 0)
    {
      totalEnergy += hitEn;
      tDistance = hitZ * fCellSizeZ;
      rDistance = hitRho * fCellSizeRho;
      phiDistance = hitPhi * fCellSizePhi;
      tFirstMoment += hitEn * tDistance;
      rFirstMoment += hitEn * rDistance;
      phiMean += hitEn * phiDistance;
      analysisManager->FillH1(4, tDistance, hitEn);
      analysisManager->FillH1(5, rDistance, hitEn);
      analysisManager->FillH1(10, hitType);
      if(hitEn > 0.0005)
      {  // e > 0.5 keV
        fCalEdep[numNonZeroThresholdCells] = hitEn;
        fCalRho[numNonZeroThresholdCells]  = hitRho;
        fCalPhi[numNonZeroThresholdCells]  = hitPhi;
        fCalZ[numNonZeroThresholdCells]    = hitZ;
        numNonZeroThresholdCells++;
        analysisManager->FillH1(13, std::log10(hitEn));
      }
    }
  }
  tFirstMoment /= totalEnergy;
  rFirstMoment /= totalEnergy;
  phiMean /= totalEnergy;
  analysisManager->FillH1(0, primaryEnergy / GeV);
  analysisManager->FillH1(1, totalEnergy / GeV);
  analysisManager->FillH1(2, totalEnergy / primaryEnergy);
  analysisManager->FillH1(3, fTimer.GetRealElapsed());
  analysisManager->FillH1(6, tFirstMoment);
  analysisManager->FillH1(7, rFirstMoment);
  analysisManager->FillH1(12, numNonZeroThresholdCells);
  // Resize to store only energy hits above threshold
  fCalEdep.resize(numNonZeroThresholdCells);
  fCalRho.resize(numNonZeroThresholdCells);
  fCalPhi.resize(numNonZeroThresholdCells);
  fCalZ.resize(numNonZeroThresholdCells);
  analysisManager->FillNtupleDColumn(0, primaryEnergy);
  analysisManager->FillNtupleDColumn(5, fTimer.GetRealElapsed());

  // Second loop over hits to calculate second moments
  for(size_t iHit = 0; iHit < hitsCollection->entries(); iHit++)
  {
    hit    = static_cast<Par04Hit*>(hitsCollection->GetHit(iHit));
    hitEn  = hit->GetEdep();
    hitZ   = hit->GetZid();
    hitRho = hit->GetRhoId();
    hitPhi = hit->GetPhiId();
    if(hitEn > 0)
    {
      tDistance = hitZ * fCellSizeZ;
      rDistance = hitRho * fCellSizeRho;
      phiDistance = hitPhi * fCellSizePhi;
      tSecondMoment += hitEn * std::pow(tDistance - tFirstMoment, 2);
      rSecondMoment += hitEn * std::pow(rDistance - rFirstMoment, 2);
      analysisManager->FillH1(11, phiDistance - phiMean, hitEn);
    }
  }
  tSecondMoment /= totalEnergy;
  rSecondMoment /= totalEnergy;
  analysisManager->FillH1(8, tSecondMoment);
  analysisManager->FillH1(9, rSecondMoment);
  analysisManager->AddNtupleRow();
}
