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
/// file:SteppingAction.cc
/// brief:
#include "SteppingAction.hh"

#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "DNAGeometry.hh"
#include "DNAHit.hh"
#include "UtilityFunctions.hh"
#include "ChromosomeMapper.hh"
#include "OctreeNode.hh"
#include "PlacementVolumeInfo.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"
#include "G4StepPoint.hh"
#include "G4VTouchable.hh"
#include "Randomize.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* evt)
  : G4UserSteppingAction()
  , fEventAction(evt)
{
  const auto* det = dynamic_cast<const DetectorConstruction*>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  fDNAGeometry = det->GetDNAGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep == 0)
    return;

  // Work out if we are in the Logical volume assigned to DNA
  const G4VTouchable* touch = aStep->GetPreStepPoint()->GetTouchable();

  // remember the chem volume is the "real" volume.
  G4LogicalVolume* dnaVol = fDNAGeometry->GetDNAChemVolumePointer();
  G4int depth             = touch->GetHistoryDepth();
  G4bool isInDNARegion    = false;

  if(depth >= 0)
  {
    fEventAction->AddCellEdep(edep);  // dousatsu
  }

  while((depth >= 0) && !isInDNARegion)
  {
    if(touch->GetVolume(depth) != nullptr)
    {
      isInDNARegion = isInDNARegion ||
                      (touch->GetVolume(depth)->GetLogicalVolume() == dnaVol);
    }
    else
    {
      G4cout << "Null pointer volume in mother with depth of "
             << touch->GetHistoryDepth() << " looking for volume at depth "
             << depth << G4endl;
      assert(depth < 0);
    }
    depth--;
  }

  if(isInDNARegion)
  {
    ProcessDNARegionHit(aStep);
  }

  G4ThreeVector prepos = aStep->GetPreStepPoint()->GetPosition();
  if(fDNAGeometry->GetChromosomeMapper()->GetChromosome(prepos) != nullptr)
  {
    DoChromosomeRegionHit(prepos, edep, false);
  }

  if(aStep->GetTrack()->GetTrackID() == 1)
  {
    fEventAction->SetPrimStopPos(prepos);
    G4double tl = aStep->GetStepLength();
    fEventAction->AddTrackLengthCell(tl);
    if(fDNAGeometry->GetChromosomeMapper()->GetChromosome(prepos) != nullptr)
    {
      fEventAction->AddTrackLengthChro(tl);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::ProcessDNARegionHit(const G4Step* aStep)
{
  G4double edep             = aStep->GetTotalEnergyDeposit();
  const G4VTouchable* touch = aStep->GetPreStepPoint()->GetTouchable();
  G4ThreeVector pos =
    aStep->GetPreStepPoint()->GetPosition() +
    G4UniformRand() * (aStep->GetPostStepPoint()->GetPosition() -
                       aStep->GetPreStepPoint()->GetPosition());
  G4VPhysicalVolume* pv = aStep->GetPreStepPoint()->GetPhysicalVolume();
  G4LogicalVolume* lv   = pv->GetLogicalVolume();

  std::vector<G4VPhysicalVolume*> molvec;
  OctreeNode* parentNode;
  G4AffineTransform transform = touch->GetHistory()->GetTopTransform();
  parentNode                  = fDNAGeometry->GetTopOctreeNode(lv);
  if(parentNode == nullptr)
  {  // See if parent has an octree.
    parentNode = fDNAGeometry->GetTopOctreeNode(pv->GetMotherLogical());
    if(parentNode != nullptr)
    {  // If mother is a node, update it to be the PV also
      // These methods are unsafe if the pv is the top volume so they
      // need to be in this if statement.
      // lv = pv->GetMotherLogical();
      pv = touch->GetVolume(1);
      transform =
        touch->GetHistory()->GetTransform(touch->GetHistoryDepth() - 1);
    }
  }

  G4ThreeVector localPos = transform.TransformPoint(pos);

  if(parentNode != nullptr)
  {  // Get Molecule Vector
    molvec = parentNode->SearchOctree(
      localPos, fDNAGeometry->GetDirectInteractionRange());
  }

  // debug output for finding logical volumes.
  if(fDNAGeometry->GetVerbosity() > 1)
  {
    if(parentNode != nullptr)
    {
      G4String lvname = pv->GetLogicalVolume()->GetName();
      G4cout << "Found octree for logical volume: " << lvname << G4endl;
    }
    else
    {
      // Don't bother doing error printouts if the particle
      // is at the DNA mother volume level.
      // Remember the Chemistry Volume is the "real" world
      if(pv->GetLogicalVolume() != fDNAGeometry->GetDNAChemVolumePointer())
      {
        G4String lvname = pv->GetLogicalVolume()->GetName();
        G4String motherlvname;
        if(lv)
          motherlvname = lv->GetName();
        G4cout << "Could not find octree for logical volume." << G4endl
               << "Particle in LV: " << lvname
               << ", with mother LV: " << motherlvname << G4endl;
      }
    }
  }

  if(!molvec.empty())
  {
    // NOTE: This loop can be optimised to avoid checking pairs twice
    G4VPhysicalVolume* closest = molvec[0];
    G4double dist              = (closest->GetTranslation() - localPos).mag();
    if(molvec.size() > 1)
    {
      G4double newdist;
      for(auto it = molvec.begin() + 1; it != molvec.end(); it++)
      {
        newdist = ((*it)->GetTranslation() - localPos).mag();
        if(newdist < dist)
        {
          dist    = newdist;
          closest = (*it);
        }
      }
    }
    G4bool isDNAHit = (dist < fDNAGeometry->GetDirectInteractionRange());
    if(isDNAHit)
    {
      DoChromosomeDNAHit(pos, localPos, edep, dist, closest->GetName(),
                         pv->GetName());
      DoChromosomeRegionHit(pos, edep, true);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::DoChromosomeRegionHit(const G4ThreeVector& pos,
                                           const G4double& edep,
                                           const G4bool& dnahit)
{
  // TODO: chromosome should probably be based on position of center of
  //       parent volume, not of current pos
  G4String chromosome =
    fDNAGeometry->GetChromosomeMapper()->GetCurrentChromosomeKey(pos);
  if(!(chromosome.empty()))
  {
    fEventAction->AddChromosomeEdep(chromosome, edep, dnahit);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::DoChromosomeDNAHit(
  const G4ThreeVector& pos, const G4ThreeVector& localPos, const G4double& edep,
  const G4double& dist, const G4String& pvName, const G4String& motherName)
{
  G4int placementIdx = std::stoi(utility::Split(motherName, '-').at(1));
  std::vector<G4String> pvNameVec = utility::Split(pvName, '-');
  molecule mol                    = utility::GetMoleculeEnum(pvNameVec.at(0));
  G4int chainIdx                  = std::stoi(pvNameVec.at(1));
  G4int strandIdx                 = std::stoi(pvNameVec.at(2));
  int64_t baseIdx                 = std::stoll(pvNameVec.at(3));

  // Convert to local chain and base index
  chainIdx = fDNAGeometry->GetGlobalChain(placementIdx, chainIdx);
  baseIdx += fDNAGeometry->GetStartIdx(placementIdx, chainIdx);

  if(fDNAGeometry->GetStrandsFlipped(placementIdx))
  {
    strandIdx = (strandIdx + 1) % 2;
  }
  G4String chromosome =
    fDNAGeometry->GetChromosomeMapper()->GetCurrentChromosomeKey(pos);
  const DNAHit* dnaHit =
    new DNAHit(mol, placementIdx, chainIdx, strandIdx, baseIdx, pos, localPos,
               edep, dist, chromosome);

  fEventAction->AddDNAHit(dnaHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

