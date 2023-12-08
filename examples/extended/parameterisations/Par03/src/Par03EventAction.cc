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
#include "Par03EventAction.hh"
#include "Par03Hit.hh"
#include "Par03DetectorConstruction.hh"

#include "G4AnalysisManager.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"

Par03EventAction::Par03EventAction(Par03DetectorConstruction* aDetector)
  : G4UserEventAction()
  , fHitCollectionID(-1)
  , fTimer()
  , fDetector(aDetector)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par03EventAction::~Par03EventAction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par03EventAction::BeginOfEventAction(const G4Event*) { fTimer.Start(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par03EventAction::EndOfEventAction(const G4Event* aEvent)
{
  fTimer.Stop();
  // Get hits collection ID (only once)
  if(fHitCollectionID == -1)
  {
    fHitCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("hits");
  }
  // Get hits collection
  auto hitsCollection = static_cast<Par03HitsCollection*>(
    aEvent->GetHCofThisEvent()->GetHC(fHitCollectionID));

  if(hitsCollection == nullptr)
  {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << fHitCollectionID;
    G4Exception("Par03EventAction::GetHitsCollection()", "MyCode0001",
                FatalException, msg);
  }
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  // Retrieve only once detector dimensions
  if(fCellSizeZ == 0)
  {
    fCellSizeZ   = fDetector->GetLength() / fDetector->GetNbOfLayers();
    fCellSizeRho = fDetector->GetRadius() / fDetector->GetNbOfRhoCells();
  }

  // Retrieve information from primary vertex and primary particle
  // To calculate shower axis and entry point to the detector
  auto primaryVertex = G4EventManager::GetEventManager()
                         ->GetConstCurrentEvent()
                         ->GetPrimaryVertex();
  auto primaryParticle   = primaryVertex->GetPrimary(0);
  G4double primaryEnergy = primaryParticle->GetTotalEnergy();
  // Estimate from vertex and particle direction the entry point to the detector
  // Calculate entrance point to the detector located at z = 0
  auto primaryDirection = primaryParticle->GetMomentumDirection();
  auto primaryEntrance  = primaryVertex->GetPosition() -
                         primaryVertex->GetPosition().z() * primaryDirection;
  G4double cosDirection = std::cos(primaryDirection.theta());
  G4double sinDirection = std::sin(primaryDirection.theta());

  // Fill histograms
  Par03Hit* hit        = nullptr;
  G4double hitEn       = 0;
  G4double totalEnergy = 0;
  G4int hitZ           = -1;
  G4int hitRho         = -1;
  G4int hitType        = -1;
  G4double tDistance = 0., rDistance = 0.;
  G4double tFirstMoment = 0., tSecondMoment = 0.;
  G4double rFirstMoment = 0., rSecondMoment = 0.;
  for(size_t iHit = 0; iHit < hitsCollection->entries(); iHit++)
  {
    hit     = static_cast<Par03Hit*>(hitsCollection->GetHit(iHit));
    hitZ    = hit->GetZid();
    hitRho  = hit->GetRhoId();
    hitEn   = hit->GetEdep();
    hitType = hit->GetType();
    if(hitEn > 0)
    {
      totalEnergy += hitEn;
      tDistance =
        hitZ * fCellSizeZ * cosDirection +
        (hitRho * fCellSizeRho - primaryEntrance.perp()) * sinDirection;
      rDistance =
        hitZ * fCellSizeZ * (-sinDirection) +
        (hitRho * fCellSizeRho - primaryEntrance.perp()) * cosDirection;
      tFirstMoment += hitEn * tDistance;
      rFirstMoment += hitEn * rDistance;
      analysisManager->FillH1(4, tDistance, hitEn);
      analysisManager->FillH1(5, rDistance, hitEn);
      analysisManager->FillH1(10, hitType);
    }
  }
  tFirstMoment /= totalEnergy;
  rFirstMoment /= totalEnergy;
  analysisManager->FillH1(0, primaryEnergy / GeV);
  analysisManager->FillH1(1, totalEnergy / GeV);
  analysisManager->FillH1(2, totalEnergy / primaryEnergy);
  analysisManager->FillH1(3, fTimer.GetRealElapsed());
  analysisManager->FillH1(6, tFirstMoment);
  analysisManager->FillH1(7, rFirstMoment);

  // Second loop over hits to calculate second moments
  for(size_t iHit = 0; iHit < hitsCollection->entries(); iHit++)
  {
    hit    = static_cast<Par03Hit*>(hitsCollection->GetHit(iHit));
    hitEn  = hit->GetEdep();
    hitZ   = hit->GetZid();
    hitRho = hit->GetRhoId();
    if(hitEn > 0)
    {
      tDistance = hitZ * fCellSizeZ * cosDirection +
                  (hitRho * fCellSizeRho - primaryEntrance.r()) * sinDirection;
      rDistance = hitZ * fCellSizeZ * (-sinDirection) +
                  (hitRho * fCellSizeRho - primaryEntrance.r()) * cosDirection;
      tSecondMoment += hitEn * std::pow(tDistance - tFirstMoment, 2);
      rSecondMoment += hitEn * std::pow(rDistance - rFirstMoment, 2);
    }
  }
  tSecondMoment /= totalEnergy;
  rSecondMoment /= totalEnergy;
  analysisManager->FillH1(8, tSecondMoment);
  analysisManager->FillH1(9, rSecondMoment);
}