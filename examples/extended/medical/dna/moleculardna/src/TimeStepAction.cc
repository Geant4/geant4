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
#include "TimeStepAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "DNAGeometry.hh"
#include "DNAHit.hh"
#include "ChromosomeMapper.hh"
#include "OctreeNode.hh"
#include "G4RunManager.hh"
#include "G4Molecule.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4ITTrackHolder.hh"
#include "G4ITTrackingManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TimeStepAction::TimeStepAction(EventAction* event)
  : G4UserTimeStepAction()
  , fEventAction(event)
  , fRadicalKillDistance(4.5 * nm)
  , fpChemistryTrackHolder(G4ITTrackHolder::Instance())
{
  AddTimeStep(1 * picosecond, 0.5 * nanosecond);
  // ctor
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction::~TimeStepAction() = default;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::StartProcessing()
{
  auto det = dynamic_cast<const DetectorConstruction*>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  fDNAGeometry = det->GetDNAGeometry();
  if(fDNAGeometry == nullptr)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "fDNAGeometry is null";
    G4Exception("TimeStepAction"
                "StartProcessing",
                "no fDNAGeometry", FatalException, exceptionDescription);
  }

  fRadicalKillDistance = fDNAGeometry->GetRadicalKillDistance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::UserPostTimeStepAction() { RadicalKillDistance(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::UserPreTimeStepAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::UserReactionAction(const G4Track&, const G4Track&,
                                        const std::vector<G4Track*>*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::RadicalKillDistance()
{
  if(fpChemistryTrackHolder == nullptr)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "fpChemistryTrackHolder is null";
    G4Exception("TimeStepAction"
                "RadicalKillDistance",
                "NO_fpChemistryTrackHolder", FatalException,
                exceptionDescription);
  }
  G4Track* trackToKill;
  G4TrackManyList::iterator it_begin =
    fpChemistryTrackHolder->GetMainList()->begin();
  G4TrackManyList::iterator it_end =
    fpChemistryTrackHolder->GetMainList()->end();
  while(it_begin != it_end)
  {
    trackToKill = nullptr;

    const G4VTouchable* touchable = it_begin->GetTouchable();
    if(touchable == nullptr)
    {
      ++it_begin;
      continue;
    }

    G4LogicalVolume* trackLV  = touchable->GetVolume()->GetLogicalVolume();
    G4LogicalVolume* motherLV = touchable->GetVolume()->GetMotherLogical();

    OctreeNode* octree_track  = fDNAGeometry->GetTopOctreeNode(trackLV);
    OctreeNode* octree_mother = fDNAGeometry->GetTopOctreeNode(motherLV);

    if((octree_track == nullptr) && (octree_mother == nullptr))
    {
      trackToKill = *it_begin;
    }
    else if(octree_track != nullptr)
    {
      const G4AffineTransform& trans =
        it_begin->GetTouchable()->GetHistory()->GetTopTransform();
      G4ThreeVector pos = trans.TransformPoint(it_begin->GetPosition());
      size_t n = octree_track->SearchOctree(pos, fRadicalKillDistance).size();
      if(n == 0)
      {
        trackToKill = *it_begin;
      }
      if(fDNAGeometry->IsInsideHistone(trackLV, pos))
      {
        trackToKill = *it_begin;
      }
    }
    else
    {
      const G4AffineTransform& trans =
        it_begin->GetTouchable()->GetHistory()->GetTopTransform();
      G4ThreeVector pos = trans.TransformPoint(it_begin->GetPosition());
      if(fDNAGeometry->IsInsideHistone(trackLV, pos))
      {
        trackToKill = *it_begin;
      }
    }
    ++it_begin;
    if(trackToKill != nullptr)
    {
      fpChemistryTrackHolder->PushToKill(trackToKill);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
