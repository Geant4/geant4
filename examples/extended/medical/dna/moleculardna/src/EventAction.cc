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
/// file:
/// brief:
#include "EventAction.hh"
#include "AnalysisManager.hh"
#include "DetectorConstruction.hh"
#include "DNAGeometry.hh"
#include "ChromosomeMapper.hh"
#include "ChromosomeHit.hh"
#include "DNAHit.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(AnalysisManager* man)
  : G4UserEventAction()
  , fAnalysisManager(man)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::Initialize()
{
  // Get DNA Geometry
  const auto* det = dynamic_cast<const DetectorConstruction*>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  fDNAGeometry = det->GetDNAGeometry();

  // Prepare Chromo Hits Map
  std::vector<G4String> keys =
    fDNAGeometry->GetChromosomeMapper()->GetChromosomeKeys();
  for(const auto& key : keys)
  {
    uint32_t key_i = G4::hashing::crc32::Hash(key);
    fChromoHitMap[key_i] = nullptr;
  }
  fInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
  if(!fInitialized)
  {
    Initialize();
  }

  fEdepCell   = 0;
  fTraLenCell = 0;
  fTraLenChro = 0;

  // Prepare DNA Hits vector
  fDNAHits.reserve(10000);

  std::vector<G4String> keys =
      fDNAGeometry->GetChromosomeMapper()->GetChromosomeKeys();
  //ChromosomeHit* theHit;
  for(const auto& it : keys)
  {
    uint32_t key_i = G4::hashing::crc32::Hash(it);
    fChromoHitMap[key_i] = new ChromosomeHit(it);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{
  // post-analysis
  fAnalysisManager->ProcessPrimary(fprimstoppos, fTraLenCell, fTraLenChro);
  fAnalysisManager->ProcessDNAHitsVector(fDNAHits);
  fAnalysisManager->ProcessChromosomeHitMap(fChromoHitMap);
  fAnalysisManager->ProcessCellEdep(fEdepCell);

  // destruction /////////////////////////////////////////////////////////////
  for(auto& it : fChromoHitMap)
  {
    delete it.second;
  }

  for(auto& fDNAHit : fDNAHits)
  {
    delete fDNAHit;
  }
  fDNAHits.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddChromosomeEdep(const G4String& key, const G4double& e,
                                    const G4bool& isDNA)
{
  uint32_t key_i = G4::hashing::crc32::Hash(key);
  auto it = fChromoHitMap.find(key_i);
  if(it != fChromoHitMap.end())
  {
    if(isDNA)
    {
      it->second->AddDNAEdep(e);
    }
    else
    {
      it->second->AddChromosomeEdep(e);
    }
  }
  else
  {
    G4ExceptionDescription errmsg;
    errmsg << "Energy deposit in unknown chromosome: " << key << G4endl;
    G4Exception("EventAction::AddChromosomeEdep", "ERR_UNKNOWN_CHROMOSOME",
                FatalException, errmsg);

  }
}
//2600754154
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddCellEdep(const G4double& edep) { fEdepCell += edep; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddTrackLengthCell(const G4double& tl) { fTraLenCell += tl; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddTrackLengthChro(const G4double& tl) { fTraLenChro += tl; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

