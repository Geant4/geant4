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
// $Id: G4RayTrajectory.cc 104015 2017-05-08 07:28:08Z gcosmo $
//
//
//

///////////////////
//G4RayTrajectory.cc
///////////////////

#include "G4RayTrajectory.hh"
#include "G4RayTrajectoryPoint.hh"
#include "G4RayTracerSceneHandler.hh"
#include "G4Step.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4TransportationManager.hh"
#include "G4ios.hh"

G4ThreadLocal G4Allocator<G4RayTrajectory>* rayTrajectoryAllocator = 0;

G4RayTrajectory :: G4RayTrajectory()
{
  positionRecord = new std::vector<G4RayTrajectoryPoint*>;
}

G4RayTrajectory :: G4RayTrajectory(G4RayTrajectory & right)
: G4VTrajectory()
{
  positionRecord = new std::vector<G4RayTrajectoryPoint*>;
  for(size_t i=0;i<right.positionRecord->size();i++)
  {
    G4RayTrajectoryPoint* rightPoint = (G4RayTrajectoryPoint*)
				((*(right.positionRecord))[i]);
    positionRecord->push_back(new G4RayTrajectoryPoint(*rightPoint));
  }
}

G4RayTrajectory :: ~G4RayTrajectory()
{
  //positionRecord->clearAndDestroy();
  for(size_t i=0;i<positionRecord->size();i++)
  { delete (*positionRecord)[i]; }
  positionRecord->clear();
  delete positionRecord;
}

void G4RayTrajectory::AppendStep(const G4Step* aStep)
{
  G4RayTrajectoryPoint* trajectoryPoint = new G4RayTrajectoryPoint();

  trajectoryPoint->SetStepLength(aStep->GetStepLength());

  G4Navigator* theNavigator 
    = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  G4bool valid;
  G4ThreeVector theLocalNormal = theNavigator->GetLocalExitNormal(&valid);
  if(valid) { theLocalNormal = -theLocalNormal; }
  G4ThreeVector theGrobalNormal 
    = theNavigator->GetLocalToGlobalTransform().TransformAxis(theLocalNormal);
  trajectoryPoint->SetSurfaceNormal(theGrobalNormal);

  G4VisManager* visManager = G4VisManager::GetInstance();
  G4RayTracerSceneHandler* sceneHandler =
  static_cast<G4RayTracerSceneHandler*>(visManager->GetCurrentSceneHandler());
  const auto& sceneVisAttsMap = sceneHandler->GetSceneVisAttsMap();

  // Make a path from the preStepPoint touchable
  G4StepPoint* preStepPoint = aStep -> GetPreStepPoint();
  const G4VTouchable* preTouchable = preStepPoint->GetTouchable();
  G4int preDepth = preTouchable->GetHistoryDepth();
  G4ModelingParameters::PVPointerCopyNoPath localPrePVPointerCopyNoPath;
  for (G4int i = preDepth; i >= 0; --i) {
    localPrePVPointerCopyNoPath.push_back
    (G4ModelingParameters::PVPointerCopyNo
     (preTouchable->GetVolume(i),preTouchable->GetCopyNumber(i)));
  }

  // Pick up the vis atts, if any, from the scene handler
  auto preIterator = sceneVisAttsMap.find(localPrePVPointerCopyNoPath);
  const G4VisAttributes* preVisAtts;
  if (preIterator != sceneVisAttsMap.end()) {
    preVisAtts = &preIterator->second;
  } else {
    preVisAtts = 0;
  }
  trajectoryPoint->SetPreStepAtt(preVisAtts);

  // Make a path from the postStepPoint touchable
  G4StepPoint* postStepPoint = aStep -> GetPostStepPoint();
  const G4VTouchable* postTouchable = postStepPoint->GetTouchable();
  G4int postDepth = postTouchable->GetHistoryDepth();
  G4ModelingParameters::PVPointerCopyNoPath localPostPVPointerCopyNoPath;
  for (G4int i = postDepth; i >= 0; --i) {
    localPostPVPointerCopyNoPath.push_back
    (G4ModelingParameters::PVPointerCopyNo
     (postTouchable->GetVolume(i),postTouchable->GetCopyNumber(i)));
  }

  // Pick up the vis atts, if any, from the scene handler
  auto postIterator = sceneVisAttsMap.find(localPostPVPointerCopyNoPath);
  const G4VisAttributes* postVisAtts;
  if (postIterator != sceneVisAttsMap.end()) {
    postVisAtts = &postIterator->second;
  } else {
    postVisAtts = 0;
  }
  trajectoryPoint->SetPostStepAtt(postVisAtts);

  positionRecord->push_back(trajectoryPoint);
}

void G4RayTrajectory::ShowTrajectory(std::ostream&) const
{ }

void G4RayTrajectory::MergeTrajectory(G4VTrajectory* secondTrajectory)
{
  if(!secondTrajectory) return;

  G4RayTrajectory* seco = (G4RayTrajectory*)secondTrajectory;
  G4int ent = seco->GetPointEntries();
  for(G4int i=0;i<ent;i++)
  { positionRecord->push_back((G4RayTrajectoryPoint*)seco->GetPoint(i)); }
  seco->positionRecord->clear();
}

