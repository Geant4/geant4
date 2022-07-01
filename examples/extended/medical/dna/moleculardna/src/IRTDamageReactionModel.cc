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
#include "IRTDamageReactionModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4Molecule.hh"
#include "G4ErrorFunction.hh"
#include <vector>
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "DNAGeometry.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4Scheduler.hh"
#include "G4UnitsTable.hh"
#include "EventAction.hh"
#include "G4EventManager.hh"
#include "DNAHit.hh"
#include "G4IRTUtils.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IRTDamageReactionModel::IRTDamageReactionModel(const G4String& name)
  : G4VDNAHitModel(name)
  , fMolecularReactionTable(G4DNAMolecularReactionTable::Instance())
{
  auto det = dynamic_cast<const DetectorConstruction*>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  fpDNAGeometry = det->GetDNAGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double IRTDamageReactionModel::GetTimeToEncounter(
  const G4MolecularConfiguration* molA, const G4MolecularConfiguration* molB,
  const G4double& distance)
{
  G4double irt          = -1;
  const auto pMoleculeA = molA;
  const auto pMoleculeB = molB;
  auto reactionData =
    fMolecularReactionTable->GetReactionData(pMoleculeA, pMoleculeB);
  G4double D =
    molA->GetDiffusionCoefficient() + molB->GetDiffusionCoefficient();
  auto kobs = reactionData->GetObservedReactionRateConstant();
  if(D == 0)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "D = 0"
                         << " for : " << molA->GetName() << " and "
                         << molB->GetName();
    G4Exception("IRTDamageReactionModel"
                "::GetTimeToEncounter()",
                "MolecularIRTDamageReactionModel002", FatalException,
                exceptionDescription);
    return -1;
  }
  G4double Reff = kobs / (4 * CLHEP::pi * D * CLHEP::Avogadro);

  if(distance < Reff)
  {
    return 0;  //
  }

  G4double Winf = 0;

  if(distance != 0)
  {
    Winf = Reff / distance;
  }
  else
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "distance = " << distance << " is uncorrected with "
                         << " Reff = " << Reff << " for : " << molA->GetName()
                         << " and " << molB->GetName();
    G4Exception("IRTDamageReactionModel"
                "::GetTimeToEncounter()",
                "MolecularIRTDamageReactionModel001", FatalException,
                exceptionDescription);
  }

  G4double U = G4UniformRand();

  if(Winf != 0 && U < Winf)
  {
    irt = (1.0 / (4 * D)) *
          std::pow((distance - Reff) / G4ErrorFunction::erfcInv(U / Winf), 2);
  }
  return irt;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IRTDamageReactionModel::RecordDNADamage()
  const  // let's suppose that this is not so much
{
  if(fpDNAPhyVolume == nullptr || fpTrack == nullptr)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "fpDNAPhyVolume == nullptr or fpTrack == nullptr";
    G4Exception("IRTDamageReactionModel"
                "RecordDNA",
                "NO_VOLUME001", FatalException, exceptionDescription);
  }
  const G4VTouchable* touchable = fpTrack->GetTouchable();
  if(touchable == nullptr)
  {
    return;
  }
  auto pPhyVolum = const_cast<G4VPhysicalVolume*>(fpDNAPhyVolume);
  const G4MolecularConfiguration* radical =
    GetMolecule(fpTrack)->GetMolecularConfiguration();

  // particle position
  const G4ThreeVector& pos_particle = fpTrack->GetPosition();
  G4AffineTransform transform     = touchable->GetHistory()->GetTopTransform();
  G4ThreeVector localpos_particle = transform.TransformPoint(pos_particle);

  // DNA position
  G4ThreeVector localPos_DNA = pPhyVolum->GetTranslation();
  G4ThreeVector globalPos_DNA =
    touchable->GetHistory()->GetTopTransform().Inverse().TransformPoint(
      localPos_DNA);

  const int64_t idx = fpDNAGeometry->GetGlobalUniqueID(pPhyVolum, touchable);

  int64_t bp         = fpDNAGeometry->GetBasePairFromUniqueID(idx);
  G4int chainID      = fpDNAGeometry->GetChainIDFromUniqueID(idx);
  G4int strandID     = fpDNAGeometry->GetStrandIDFromUniqueID(idx);
  molecule mol       = fpDNAGeometry->GetMoleculeFromUniqueID(idx);
  G4int placementIdx = fpDNAGeometry->GetPlacementIndexFromUniqueID(idx);

  G4String chromosome =
    fpDNAGeometry->GetChromosomeMapper()->GetCurrentChromosomeKey(
      globalPos_DNA);

  auto dnaHit = new DNAHit(mol, placementIdx, chainID, strandID, bp,
                           globalPos_DNA, localPos_DNA, chromosome, radical);

  dynamic_cast<EventAction*>(
    G4EventManager::GetEventManager()->GetUserEventAction())
    ->AddDNAHit(dnaHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IRTDamageReactionModel::MakeReaction(const G4Track& track)
{
  // G4Track *pTrackA = const_cast<G4Track *>(&track);
  if(track.GetTrackStatus() == fStopAndKill)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "this track is killed";
    G4Exception("IRTDamageReactionModel"
                "MakeReaction",
                "NO_TRACK02", FatalException, exceptionDescription);
  }
  if(G4Scheduler::Instance()->GetVerbose() != 0)
  {
    G4cout << "At time : " << std::setw(7)
           << G4BestUnit(G4Scheduler::Instance()->GetGlobalTime(), "Time")
           << " Reaction : " << GetIT(track)->GetName() << " ("
           << track.GetTrackID() << ") + " << fpDNAPhyVolume->GetName()
           << " -> ";
    G4cout << " Damaged " + fpDNAPhyVolume->GetName();
    G4cout << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// DoReaction means : kill species and record DNA damage
G4bool IRTDamageReactionModel::DoReaction(const G4Track& track,
                                          const G4double& reactionTime,
                                          const DNANode& vp)
{
  fReactionTime = reactionTime;

  if(fReactionTime == G4Scheduler::Instance()->GetLimitingTimeStep())
  {
    return false;
  }

  fpTrack        = &track;
  fpDNAPhyVolume = std::get<const G4VPhysicalVolume*>(vp);
  MakeReaction(track);
  RecordDNADamage();
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double IRTDamageReactionModel::CalculateReactionTime(const G4Track& track,
                                                       DNANode& vp)
{
  fpTrack        = nullptr;
  fminTimeStep   = DBL_MAX;
  fReactionTime  = DBL_MAX;
  fpDNAPhyVolume = nullptr;
  if(fpDNAGeometry == nullptr)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "no fpDNAGeometry" << G4endl;
    G4Exception("IRTDamageReactionModel"
                "::CalculateReactionTime()",
                "MolecularIRTDamageReactionModel007", FatalException,
                exceptionDescription);
  }

  auto pMoleculeA          = GetMolecule(track);
  auto pMolConfA           = pMoleculeA->GetMolecularConfiguration();
  const auto pReactantList = fMolecularReactionTable->CanReactWith(pMolConfA);

  if(pReactantList == nullptr)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "!!!!!!!!!!!!!!!!!!!!" << G4endl;
    G4cout << "!!! WARNING" << G4endl;
    G4cout
      << "IRTDamageReactionModel::CalculateReactionTime will return infinity "
         "for the reaction because the molecule "
      << pMoleculeA->GetName()
      << " does not have any reactants given in the reaction table." << G4endl;
    G4cout << "!!!!!!!!!!!!!!!!!!!!" << G4endl;
    G4Exception("IRTDamageReactionModel"
                "::CalculateReactionTime()",
                "MolecularIRTDamageReactionModel003", FatalException,
                exceptionDescription);
    return -1;
  }

  size_t nbReactives = pReactantList->size();

  if(nbReactives == 0)
  {
    G4cout << "!!!!!!!!!!!!!!!!!!!!" << G4endl;
    G4cout << "!!! WARNING" << G4endl;
    G4cout
      << "IRTDamageReactionModel::CalculateReactionTime will return infinity "
         "for the reaction because the molecule "
      << pMoleculeA->GetName()
      << " does not have any reactants given in the reaction table."
      << "This message can also result from a wrong implementation of the "
         "reaction table."
      << G4endl;
    G4cout << "!!!!!!!!!!!!!!!!!!!!" << G4endl;
    return -1;
  }
  const G4VTouchable* touchable = track.GetTouchable();
  if(touchable == nullptr)
  {
    return -1;
  }

  const G4LogicalVolume* logicalVolume =
    touchable->GetVolume()->GetLogicalVolume();
  const G4ThreeVector& globalPos_track = track.GetPosition();
  const G4ThreeVector& localPos =
    touchable->GetHistory()->GetTopTransform().TransformPoint(globalPos_track);

  G4double D = pMoleculeA->GetDiffusionCoefficient();
  std::vector<G4VPhysicalVolume*> result_pv;
  result_pv.clear();
  fpDNAGeometry->FindNearbyMolecules(logicalVolume, localPos, result_pv,
                                     G4IRTUtils::GetRCutOff(100 * ns));

  if(result_pv.empty())
  {
    return -1;
  }
  for(const auto& physicalVolume : result_pv)
  {
    const G4Material* material =
      physicalVolume->GetLogicalVolume()->GetMaterial();

    G4MolecularConfiguration* dna_molConf =
      G4DNAMolecularMaterial::Instance()->GetMolecularConfiguration(material);
    auto it =
      std::find(pReactantList->begin(), pReactantList->end(), dna_molConf);
    if(it == pReactantList->end())
    {
      continue;
    }

    G4ThreeVector localPos_DNA = physicalVolume->GetTranslation();
    G4ThreeVector globalPos_DNA =
      touchable->GetHistory()->GetTopTransform().Inverse().TransformPoint(
        localPos_DNA);
    G4double distance = (localPos_DNA - localPos).mag();

    G4double distance2 = distance * distance;
    auto reactionData =
      G4DNAMolecularReactionTable::Instance()->GetReactionData(pMolConfA,
                                                               dna_molConf);

    G4double kobs = reactionData->GetObservedReactionRateConstant();
    G4double Reff = kobs / (4 * CLHEP::pi * D * CLHEP::Avogadro);

    if(distance2 < Reff * Reff)
    {
      fminTimeStep = 0.;
      vp           = physicalVolume;
    }
    else
    {
      G4double tempMinET = GetTimeToEncounter(pMolConfA, dna_molConf, distance);
      if(tempMinET > G4Scheduler::Instance()->GetEndTime() || tempMinET < 0)
      {
        continue;
      }
      if(tempMinET >= fminTimeStep)
      {
        continue;
      }
      fminTimeStep = tempMinET;
      vp           = physicalVolume;
    }
  }
  if(fminTimeStep > G4Scheduler::Instance()->GetLimitingTimeStep() &&
     fminTimeStep < G4Scheduler::Instance()->GetEndTime())
  {
    fminTimeStep = G4Scheduler::Instance()->GetLimitingTimeStep();
  }
  return fminTimeStep;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
