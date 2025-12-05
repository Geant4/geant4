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
/// \file TrackerSD.cc
/// \brief Implementation of the TrackerSD class

#include "TrackerSD.hh"

#include "G4AnalysisManager.hh"
#include "G4Exception.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessType.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4VProcess.hh"
#include "Randomize.hh"
#include "globals.hh"

#include "DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TrackerSD::TrackerSD(const G4String& sdName, const G4String& hitsCollectionName,
                     const int /*depthIndex*/)
  : G4VSensitiveDetector(sdName)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::~TrackerSD() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection =
    new TrackerHitColl(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  G4int collID =
    G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);

  hce->AddHitsCollection(collID, fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  // Get the track
  G4Track* aTrack = aStep->GetTrack();

  // Get Pre and Post step points
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint* postStepPoint = aStep->GetPostStepPoint();

  // Get ParentID
  G4int parentID = aTrack->GetParentID();

  // Only for primary particles
  if (parentID == 0) {
    // Check if the particle is entering or exiting
    if (preStepPoint->GetStepStatus() == fGeomBoundary) {
      // Get the pre-step energy
      G4double preEnergy = preStepPoint->GetKineticEnergy();

      // Fill the histogram
      G4AnalysisManager::Instance()->FillH1(12, preEnergy);
    }

    if (postStepPoint->GetStepStatus() == fGeomBoundary) {
      // Get the post-step energy
      G4double postEnergy = postStepPoint->GetKineticEnergy();

      // Fill the histogram
      G4AnalysisManager::Instance()->FillH1(13, postEnergy);
    }
  }

  // Get the energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();
  if (edep == 0.) return false;

  // Declare a new hit
  auto newHit = new TrackerHit();

  newHit->SetEdep(edep);

  // Get the position of the hit
  G4ThreeVector Pos = postStepPoint->GetPosition();
  newHit->SetPosition(Pos);

  // Add the hit to the collection
  fHitsCollection->insert(newHit);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  const auto nofHits = fHitsCollection->entries();

  if (verboseLevel > 0) {
    G4cout << "  Hits Collection: in this event, we have " << nofHits
           << " hits in the tracker." << G4endl;
  }

  // Get the DetectorConstruction instance
  const DetectorConstruction* detConstruction =
    static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  // Declare variables for the dimensions of the detector
  G4double HitSelRegZ = 0.;
  G4double HitSelRegXY = 0.;
  G4double site_radius = 0.;
  G4double DetDensity = 0.;
  G4ThreeVector hitPos;
  G4ThreeVector randCenterPos;

  if (detConstruction) {
    site_radius = detConstruction->GetSiteRadius();
    HitSelRegZ = detConstruction->GetHitSelRegZ();
    HitSelRegXY = detConstruction->GetHitSelRegXY();
    auto mat = G4Material::GetMaterial(detConstruction->GetMaterial());
    if (mat)
      DetDensity = mat->GetDensity();
    else {
      G4ExceptionDescription msg;
      msg << "Material not found."
          << "Something unexpected has occurred." << G4endl;
      G4Exception("TrackerSD::EndOfEvent()", "TrackerSD002", JustWarning, msg);
    }
  }
  else {
    G4ExceptionDescription msg;
    msg << "Detector construction not found."
        << "Something unexpected has occurred." << G4endl;
    G4Exception("TrackerSD::EndOfEvent()", "TrackerSD001", JustWarning, msg);
  }
  G4int nHsel = 0;  // number of hits in the selection region
  G4int nHsite = 0;  // number of hits in the site
  G4int nHint = 0;  // number of hits both in site and selection region
  G4double evtEdep = 0.;  // energy deposit in the event

  if (nofHits > 0) {
    const G4int maxTries = 1000;
    G4int tries = 0;
    G4bool found = false;

    while (tries < maxTries && !found) {
      std::size_t randHit = static_cast<std::size_t>(G4UniformRand() * nofHits);
      auto hit = (*fHitsCollection)[randHit];
      G4ThreeVector hitPosition = hit->GetPosition();

      if (hitPosition.x() > -HitSelRegXY / 2
          && hitPosition.x() < HitSelRegXY / 2
          && hitPosition.y() > -HitSelRegXY / 2
          && hitPosition.y() < HitSelRegXY / 2
          && hitPosition.z() > -HitSelRegZ / 2
          && hitPosition.z() < HitSelRegZ / 2)
      {
        hitPos = hitPosition;  // store valid hit position

        found = true;
      }
      ++tries;
    }

    if (found) {
      G4double site_radius2 = site_radius * site_radius;

      // Random placement of a sphere center containing the hit
      G4double xRand, yRand, zRand, randRad2;
      do {
        xRand = (2 * G4UniformRand() - 1) * site_radius;
        yRand = (2 * G4UniformRand() - 1) * site_radius;
        zRand = (2 * G4UniformRand() - 1) * site_radius;
        randRad2 = xRand * xRand + yRand * yRand + zRand * zRand;
      } while (randRad2 > site_radius2);

      randCenterPos = G4ThreeVector(xRand + hitPos.x(), yRand + hitPos.y(),
                                    zRand + hitPos.z());
    }

    else {
      G4cout
        << "In this event, no hits were found in the hit selection "
        << "region after trying " << tries << " times." << G4endl
        << "Please consider increasing the size of the hit selection region."
        << G4endl;
      return;  // Skip analysis if no valid hit was found
    }
  }

  // Loop over hits
  for (std::size_t jj = 0; jj < nofHits; jj++) {
    auto hit2 = (*fHitsCollection)[jj];
    G4ThreeVector hit2Position = hit2->GetPosition();

    G4bool inSlab =
      (hit2Position.x() > -HitSelRegXY / 2 && hit2Position.x() < HitSelRegXY / 2
       && hit2Position.y() > -HitSelRegXY / 2
       && hit2Position.y() < HitSelRegXY / 2
       && hit2Position.z() > -HitSelRegZ / 2
       && hit2Position.z() < HitSelRegZ / 2);

    // Check if the hit is within the site
    G4double dist = (hit2Position - randCenterPos).mag();
    if (dist < site_radius) {
      nHsite++;
      evtEdep += hit2->GetEdep();
    }

    // Check if the hit is within the selection region
    if (inSlab) nHsel++;

    // Check if the hit is within the selection region and the site
    if (dist < site_radius && inSlab) nHint++;
  }

  // Access the analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Calculate lineal energy
  G4double y = (evtEdep) / ((4. / 3.) * site_radius);  // in keV/um

  // Calculate specific energy
  G4double mass = ((4. / 3.) * CLHEP::pi * site_radius * site_radius
                   * site_radius * DetDensity);
  G4double z = (evtEdep / mass);

  // Fill histograms
  if (nHint > 0) {
    // Define the weight
    G4double wght = G4double(nHsel) / G4double(nHint);

    // Histogram 0: Single-event energy imparted in keV
    analysisManager->FillH1(0, evtEdep / keV, wght);

    // Histogram 1: Weighted single-event energy imparted in keV
    analysisManager->FillH1(1, evtEdep / keV, (evtEdep * wght / keV));

    // Histogram 2: Squared weighted energy imparted per event in keV^2
    analysisManager->FillH1(2, evtEdep / keV,
                            ((evtEdep * evtEdep * wght) / (keV * keV)));

    // Histogram 3: Lineal energy in keV/um
    analysisManager->FillH1(3, y / (keV / um), wght);

    // Histogram 4: Dose-weighted lineal energy in keV/um
    analysisManager->FillH1(4, y / (keV / um), (y * wght / (keV / um)));

    // Histogram 5: Squared weighted lineal energy in (keV/um)^2
    analysisManager->FillH1(5, y / (keV / um),
                            (y * y * wght / ((keV / um) * (keV / um))));

    // Histogram 6: Specific energy in keV/g
    analysisManager->FillH1(6, z / gray, wght);

    // Histogram 7: Weighted specific energy in keV/g
    analysisManager->FillH1(7, z / gray, (z * wght / gray));

    // Histogram 8: Squared-weighted specific energy in (keV/g)^2
    analysisManager->FillH1(8, z / gray, (z * z * wght / (gray * gray)));

    // Histogram 9: Number of hits in the selection region
    analysisManager->FillH1(9, nHsel, 1.);

    // Histogram 10: Number of hits in the site
    analysisManager->FillH1(10, nHsite, 1.);

    // Histogram 11: Number of hits both in site and selection region
    analysisManager->FillH1(11, nHint, 1.);

    // Histogram 14: Number of hits in site vs Energy imparted per event
    analysisManager->FillH2(0, evtEdep / keV, nHsite, 1.);
  }
  else {
    G4ExceptionDescription msg;
    msg << "In this event, we had nHint == 0. "
        << "Something unexpected has occurred." << G4endl;
    G4Exception("TrackerSD::EndOfEvent()", "TrackerSD003", JustWarning, msg);
  }
}