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
// $Id: G4TrajectoryOriginVolumeFilter.cc 95115 2016-01-26 16:33:29Z gcosmo $
//
// Filter trajectories according to volume name. Only registered 
// volumes will pass the filter.
//
// Jane Tinslay May 2006
//
#include "G4TrajectoryOriginVolumeFilter.hh"
#include "G4TransportationManager.hh"
#include "G4VTrajectoryPoint.hh"

G4TrajectoryOriginVolumeFilter::G4TrajectoryOriginVolumeFilter(const G4String& name)
  :G4SmartFilter<G4VTrajectory>(name)
{}

G4TrajectoryOriginVolumeFilter::~G4TrajectoryOriginVolumeFilter() {}

bool
G4TrajectoryOriginVolumeFilter::Evaluate(const G4VTrajectory& traj) const
{
  G4VTrajectoryPoint* aTrajectoryPoint = traj.GetPoint(0);
  assert (0 != aTrajectoryPoint);

  G4Navigator* navigator =
  G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  
  G4VPhysicalVolume* volume = navigator->LocateGlobalPointAndSetup
  (aTrajectoryPoint->GetPosition(), nullptr,false,true);

  // Logical volume
  G4LogicalVolume* logicalVolume = volume->GetLogicalVolume();
  assert (0 != logicalVolume);

  // Get volume names
  G4String logicalName = logicalVolume->GetName();
  G4String physicalName = volume->GetName();

  if (GetVerbose()) {
    G4cout<<"G4TrajectoryOriginVolumeFilter processing trajectory with originating volume "<<G4endl;
    G4cout<<"logical and physical names:  "<<logicalName<<" "<<physicalName<<G4endl;
  }
  // Search for logical volume name
  std::vector<G4String>::const_iterator iterLogical = std::find(fVolumes.begin(), fVolumes.end(), logicalName);

  // Keep if logical volume registered
  if (iterLogical != fVolumes.end()) return true;

  // Repeat for physical volume name
  std::vector<G4String>::const_iterator iterPhysical = std::find(fVolumes.begin(), fVolumes.end(), physicalName);

  if (iterPhysical != fVolumes.end()) return true;

  // Volume names not registered
  return false;
}

void
G4TrajectoryOriginVolumeFilter::Add(const G4String& volume)
{
  fVolumes.push_back(volume);
}

void
G4TrajectoryOriginVolumeFilter::Print(std::ostream& ostr) const
{
  ostr<<"Volume names registered: "<<G4endl;
  std::vector<G4String>::const_iterator iter = fVolumes.begin();
  
  while (iter != fVolumes.end()) {
    ostr<<*iter<<G4endl;    
    iter++;
  }
}

void 
G4TrajectoryOriginVolumeFilter::Clear()
{
  // Clear volume vector
  fVolumes.clear();
}
