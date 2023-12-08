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

#include "Par04DetectorConstruction.hh"  // for Par04DetectorConstruction
#include "Par04Hit.hh"                   // for Par04Hit, Par04HitsCollection
#include "Par04ParallelFullWorld.hh"

#include "G4AnalysisManager.hh"          // for G4AnalysisManager
#include "G4Event.hh"                    // for G4Event
#include "G4EventManager.hh"             // for G4EventManager
#include "G4HCofThisEvent.hh"            // for G4HCofThisEvent
#include "G4SDManager.hh"                // for G4SDManager

#include <CLHEP/Units/SystemOfUnits.h>   // for GeV
#include <CLHEP/Vector/ThreeVector.h>    // for Hep3Vector
#include "G4Exception.hh"                // for G4Exception, G4ExceptionDesc...
#include "G4ExceptionSeverity.hh"        // for FatalException
#include "G4GenericAnalysisManager.hh"   // for G4GenericAnalysisManager
#include "G4PrimaryParticle.hh"          // for G4PrimaryParticle
#include "G4PrimaryVertex.hh"            // for G4PrimaryVertex
#include "G4SystemOfUnits.hh"            // for GeV
#include "G4THitsCollection.hh"          // for G4THitsCollection
#include "G4ThreeVector.hh"              // for G4ThreeVector
#include "G4Timer.hh"                    // for G4Timer
#include "G4UserEventAction.hh"          // for G4UserEventAction
#include <algorithm>                     // for max
#include <cmath>                         // for log10
#include <cstddef>                      // for size_t
#include <ostream>                       // for basic_ostream::operator<<

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04EventAction::Par04EventAction(Par04DetectorConstruction* aDetector,
                                   Par04ParallelFullWorld* aParallel)
  : G4UserEventAction()
  , fHitCollectionID(-1)
  , fPhysicalFullHitCollectionID(-1)
  , fPhysicalFastHitCollectionID(-1)
  , fTimer()
  , fDetector(aDetector)
  , fParallel(aParallel)
{
  fCellNbRho = aDetector->GetMeshNbOfCells().x();
  fCellNbPhi = aDetector->GetMeshNbOfCells().y();
  fCellNbZ   = aDetector->GetMeshNbOfCells().z();
  fCalEdep.reserve(fCellNbRho * fCellNbPhi * fCellNbZ);
  fCalRho.reserve(fCellNbRho * fCellNbPhi * fCellNbZ);
  fCalPhi.reserve(fCellNbRho * fCellNbPhi * fCellNbZ);
  fCalZ.reserve(fCellNbRho * fCellNbPhi * fCellNbZ);
  fCalPhysicalEdep.reserve(fPhysicalNbLayers * fPhysicalNbRows * fPhysicalNbSlices);
  fCalPhysicalLayer.reserve(fPhysicalNbLayers * fPhysicalNbRows * fPhysicalNbSlices);
  fCalPhysicalSlice.reserve(fPhysicalNbLayers * fPhysicalNbRows * fPhysicalNbSlices);
  fCalPhysicalRow.reserve(fPhysicalNbLayers * fPhysicalNbRows * fPhysicalNbSlices);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04EventAction::~Par04EventAction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04EventAction::BeginOfEventAction(const G4Event*) {
  StartTimer();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04EventAction::StartTimer() {
  fTimer.Start();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04EventAction::StopTimer() {
  fTimer.Stop();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04EventAction::EndOfEventAction(const G4Event* aEvent)
{
  G4SDManager::GetSDMpointer()->GetHCtable();
  StopTimer();

  // Get hits collection ID (only once)
  if(fHitCollectionID == -1)
  {
    fHitCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("hits");
  }
  if(fPhysicalFullHitCollectionID == -1)
  {
    fPhysicalFullHitCollectionID =
      G4SDManager::GetSDMpointer()->GetCollectionID("physicalCellsFullSim");
  }
  if(fPhysicalFastHitCollectionID == -1)
  {
    fPhysicalFastHitCollectionID =
      G4SDManager::GetSDMpointer()->GetCollectionID("physicalCellsFastSim");
  }
  // Get hits collection
  auto hitsCollection =
    static_cast<Par04HitsCollection*>(aEvent->GetHCofThisEvent()->GetHC(fHitCollectionID));
  auto physicalFullHitsCollection =
    static_cast<Par04HitsCollection*>(aEvent->GetHCofThisEvent()
                                      ->GetHC(fPhysicalFullHitCollectionID));
  auto physicalFastHitsCollection =
    static_cast<Par04HitsCollection*>(aEvent->GetHCofThisEvent()
                                      ->GetHC(fPhysicalFastHitCollectionID));

  if(hitsCollection == nullptr)
  {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << fHitCollectionID;
    G4Exception("Par04EventAction::GetHitsCollection()", "MyCode0001", FatalException, msg);
  }
  if(physicalFullHitsCollection == nullptr)
  {
    G4ExceptionDescription msg;
    msg << "Cannot access physical full sim hitsCollection ID " << fPhysicalFullHitCollectionID;
    G4Exception("Par04EventAction::GetHitsCollection()", "MyCode0001", FatalException, msg);
  }
  if(physicalFastHitsCollection == nullptr)
  {
    G4ExceptionDescription msg;
    msg << "Cannot access physical fast sim hitsCollection ID " << fPhysicalFastHitCollectionID;
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
  if(fPhysicalNbLayers == 0)
  {
    fPhysicalNbLayers = fParallel->GetNbOfLayers();
    fPhysicalNbSlices = fParallel->GetNbOfSlices();
    fPhysicalNbRows = fParallel->GetNbOfRows();
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
  fCalEdep.resize(fCellNbRho * fCellNbPhi * fCellNbZ,0);
  fCalRho.resize(fCellNbRho * fCellNbPhi * fCellNbZ,0);
  fCalPhi.resize(fCellNbRho * fCellNbPhi * fCellNbZ,0);
  fCalZ.resize(fCellNbRho * fCellNbPhi * fCellNbZ,0);
  fCalPhysicalEdep.resize(fPhysicalNbLayers * fPhysicalNbRows * fPhysicalNbSlices,0);
  fCalPhysicalLayer.resize(fPhysicalNbLayers * fPhysicalNbRows * fPhysicalNbSlices,0);
  fCalPhysicalSlice.resize(fPhysicalNbLayers * fPhysicalNbRows * fPhysicalNbSlices,0);
  fCalPhysicalRow.resize(fPhysicalNbLayers * fPhysicalNbRows * fPhysicalNbSlices,0);

  // Fill histograms
  Par04Hit* hit                  = nullptr;
  G4double hitEn                 = 0;
  G4double totalEnergy           = 0;
  G4int hitNum                   = 0;
  G4int totalNum                 = 0;
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
    hitNum  = hit->GetNdep();
    hitType = hit->GetType();
    if(hitEn > 0)
    {
      totalEnergy += hitEn;
      totalNum += hitNum;
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
        analysisManager->FillH1(15, hitNum);
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
  analysisManager->FillH1(14, totalNum);
  // Resize to store only energy hits above threshold
  fCalEdep.resize(numNonZeroThresholdCells);
  fCalRho.resize(numNonZeroThresholdCells);
  fCalPhi.resize(numNonZeroThresholdCells);
  fCalZ.resize(numNonZeroThresholdCells);
  analysisManager->FillNtupleDColumn(0, 0, primaryEnergy);
  analysisManager->FillNtupleDColumn(0, 1, fTimer.GetRealElapsed());
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

  // Fill ntuple with physical readout data
  G4double totalPhysicalEnergy   = 0;
  totalNum = 0;
  hitEn = 0;
  hitNum = 0;
  G4int hitLayer                 = -1;
  G4int hitRow                   = -1;
  G4int hitSlice                 = -1;
  numNonZeroThresholdCells = 0;
  for(size_t iHit = 0; iHit < physicalFullHitsCollection->entries(); iHit++)
  {
    hit     = static_cast<Par04Hit*>(physicalFullHitsCollection->GetHit(iHit));
    hitLayer = hit->GetRhoId();
    hitRow  = hit->GetZid();
    hitSlice  = hit->GetPhiId();
    hitEn   = hit->GetEdep();
    hitNum  = hit->GetNdep();
    if(hitEn > 0)
    {
      totalPhysicalEnergy += hitEn;
      totalNum += hitNum;
      if(hitEn > 0.0005)
      {  // e > 0.5 keV
        fCalPhysicalEdep[numNonZeroThresholdCells] = hitEn;
        fCalPhysicalLayer[numNonZeroThresholdCells]  = hitLayer;
        fCalPhysicalRow[numNonZeroThresholdCells]  = hitRow;
        fCalPhysicalSlice[numNonZeroThresholdCells]    = hitSlice;
        numNonZeroThresholdCells++;
        analysisManager->FillH1(19, std::log10(hitEn));
        analysisManager->FillH1(21, hitNum);
      }
    }
  }for(size_t iHit = 0; iHit < physicalFastHitsCollection->entries(); iHit++)
  {
    hit     = static_cast<Par04Hit*>(physicalFastHitsCollection->GetHit(iHit));
    hitLayer = hit->GetRhoId();
    hitRow  = hit->GetZid();
    hitSlice  = hit->GetPhiId();
    hitEn   = hit->GetEdep();
    hitNum  = hit->GetNdep();
    if(hitEn > 0)
    {
      totalPhysicalEnergy += hitEn;
      totalNum += hitNum;
      if(hitEn > 0.0005)
      {  // e > 0.5 keV
        fCalPhysicalEdep[numNonZeroThresholdCells] = hitEn;
        fCalPhysicalLayer[numNonZeroThresholdCells]  = hitLayer;
        fCalPhysicalRow[numNonZeroThresholdCells]  = hitRow;
        fCalPhysicalSlice[numNonZeroThresholdCells]    = hitSlice;
        numNonZeroThresholdCells++;
        analysisManager->FillH1(19, std::log10(hitEn));
        analysisManager->FillH1(21, hitNum);
      }
    }
  }
  analysisManager->FillH1(16, totalPhysicalEnergy / GeV);
  analysisManager->FillH1(17, totalPhysicalEnergy / primaryEnergy);
  analysisManager->FillH1(18, numNonZeroThresholdCells);
  analysisManager->FillH1(20, totalNum);
  fCalPhysicalEdep.resize(numNonZeroThresholdCells);
  fCalPhysicalLayer.resize(numNonZeroThresholdCells);
  fCalPhysicalSlice.resize(numNonZeroThresholdCells);
  fCalPhysicalRow.resize(numNonZeroThresholdCells);
  analysisManager->AddNtupleRow(0);
  analysisManager->AddNtupleRow(1);
  analysisManager->AddNtupleRow(2);
}
