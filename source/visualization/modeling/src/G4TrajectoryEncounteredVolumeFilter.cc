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
// $Id: G4TrajectoryEncounteredVolumeFilter.cc 95107 2016-01-26 12:41:11Z allison $
//
// Filter trajectories according to encountered volume name. Only registered 
// volumes will pass the filter.
//
// John Allison, Feb 2016
//
#include "G4TrajectoryEncounteredVolumeFilter.hh"
#include "G4TransportationManager.hh"
#include "G4VTrajectoryPoint.hh"
#include "G4AttValue.hh"
#include "G4RichTrajectory.hh"

G4TrajectoryEncounteredVolumeFilter::G4TrajectoryEncounteredVolumeFilter(const G4String& name)
  :G4SmartFilter<G4VTrajectory>(name)
{}

G4TrajectoryEncounteredVolumeFilter::~G4TrajectoryEncounteredVolumeFilter() {}

bool
G4TrajectoryEncounteredVolumeFilter::Evaluate(const G4VTrajectory& traj) const
{
  try
  {
    const G4RichTrajectory& richTrajectory = dynamic_cast<const G4RichTrajectory&>(traj);

    for (const auto& pvname: fVolumes) {
      for (G4int iPoint = 0; iPoint < richTrajectory.GetPointEntries(); iPoint++) {
        G4VTrajectoryPoint* point = richTrajectory.GetPoint(iPoint);
        if (!point) continue;
        std::vector<G4AttValue>* attValues = point->CreateAttValues();
        std::vector<G4AttValue>::const_iterator iAtt;
        for (iAtt = attValues->begin(); iAtt != attValues->end(); ++iAtt) {
          if (iAtt->GetName() == "PostVPath" &&
              iAtt->GetValue().contains(pvname)) break;
        }
        if (iAtt != attValues->end()) {  // Required value found
          return true;  // First found pvname determines selection.
        }
      }
    }
    return false;
  }

  catch (std::bad_cast)
  {
    G4ExceptionDescription ed;
    ed << "Requires G4RichTrajectory - \"/vis/scene/add/trajectories rich\"";
    G4Exception
    ("G4TrajectoryEncounteredVolumeFilter::Evaluate(const G4VTrajectory& traj)",
     "modeling0126",
     JustWarning, ed);
    return false;
  }
}

void
G4TrajectoryEncounteredVolumeFilter::Add(const G4String& volume)
{
  fVolumes.push_back(volume);
}

void
G4TrajectoryEncounteredVolumeFilter::Print(std::ostream& ostr) const
{
  ostr<<"Volume names registered: "<<G4endl;
  std::vector<G4String>::const_iterator iter = fVolumes.begin();
  
  while (iter != fVolumes.end()) {
    ostr<<*iter<<G4endl;    
    iter++;
  }
}

void 
G4TrajectoryEncounteredVolumeFilter::Clear()
{
  // Clear volume vector
  fVolumes.clear();
}
