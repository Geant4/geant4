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
// $Id: G4TrajectoryDrawByEncounteredVolume.cc 95107 2016-01-26 12:41:11Z allison $
//
// John Allison February 2016, based on
// G4TrajectoryDrawByVolume.hh  Jane Tinslay March 2006

#include "G4TrajectoryDrawByEncounteredVolume.hh"
#include "G4TrajectoryDrawerUtils.hh"
#include "G4VTrajectory.hh"
#include "G4VisTrajContext.hh"
#include "G4VTrajectoryPoint.hh"
#include "G4AttValue.hh"
#include "G4RichTrajectory.hh"

G4TrajectoryDrawByEncounteredVolume::G4TrajectoryDrawByEncounteredVolume
(const G4String& name, G4VisTrajContext* context)
  :G4VTrajectoryModel(name, context)
  ,fDefault(G4Colour::Grey())
{}

G4TrajectoryDrawByEncounteredVolume::~G4TrajectoryDrawByEncounteredVolume() {}

void
G4TrajectoryDrawByEncounteredVolume::Draw(const G4VTrajectory& traj, const G4bool& visible) const
{
  try
  {
    const G4RichTrajectory& richTrajectory = dynamic_cast<const G4RichTrajectory&>(traj);

    G4Colour colour(fDefault);
    G4String soughtPVName("none");

    for (const auto& item: fMap.GetBasicMap()) {
      soughtPVName = item.first;
      for (G4int iPoint = 0; iPoint < richTrajectory.GetPointEntries(); iPoint++) {
        G4VTrajectoryPoint* point = richTrajectory.GetPoint(iPoint);
        if (!point) continue;
        std::vector<G4AttValue>* attValues = point->CreateAttValues();
        std::vector<G4AttValue>::const_iterator iAtt;
        for (iAtt = attValues->begin(); iAtt != attValues->end(); ++iAtt) {
          if (iAtt->GetName() == "PostVPath" &&
              iAtt->GetValue().contains(soughtPVName)) break;
        }
        if (iAtt != attValues->end()) {  // Required value found
          fMap.GetColour(soughtPVName, colour);
          break;  // First found pvname determines colour.
        }
      }
    }

    G4VisTrajContext myContext(GetContext());

    myContext.SetLineColour(colour);
    myContext.SetVisible(visible);

    if (GetVerbose()) {
      G4cout
      << "G4TrajectoryDrawByEncounteredVolume drawer named " << Name()
      << ", drawing trajectory touching physical volume " << soughtPVName
      << ", with configuration:" << G4endl;
      myContext.Print(G4cout);
    }
    
    G4TrajectoryDrawerUtils::DrawLineAndPoints(richTrajectory, myContext);
    
  }

  catch (std::bad_cast)
  {
    G4ExceptionDescription ed;
    ed << "Requires G4RichTrajectory - \"/vis/scene/add/trajectories rich\"";
    G4Exception
    ("G4TrajectoryDrawByEncounteredVolume::Draw(const G4VTrajectory& traj,...",
     "modeling0125",
     JustWarning, ed);
    return;
  }
}

void
G4TrajectoryDrawByEncounteredVolume::SetDefault(const G4String& colour)
{
  G4Colour myColour;      

  // Will not modify default colour if colour key does not exist  
  if (!G4Colour::GetColour(colour, myColour)) {
    G4ExceptionDescription ed;
    ed << "G4Colour with key "<<colour<<" does not exist ";
    G4Exception
      ("G4TrajectoryDrawByEncounteredParticleID::SetDefault(const G4String& colour)",
       "modeling0123", JustWarning, ed);
    return;
  }

  SetDefault(myColour);
}

void
G4TrajectoryDrawByEncounteredVolume::SetDefault(const G4Colour& colour)
{
  fDefault = colour;
}

void
G4TrajectoryDrawByEncounteredVolume::Set(const G4String& pvname, const G4String& colour)
{
  fMap.Set(pvname, colour);
}

void
G4TrajectoryDrawByEncounteredVolume::Set(const G4String& pvname, const G4Colour& colour)
{
  fMap[pvname] = colour;
}

void
G4TrajectoryDrawByEncounteredVolume::Print(std::ostream& ostr) const
{
  ostr
  << "G4TrajectoryDrawByEncounteredVolume model "<< Name()
  << ", colour scheme: "
  << ", Default " << fDefault
  << std::endl;

  fMap.Print(ostr);

  ostr << "Default configuration:" << std::endl;
  GetContext().Print(ostr);
}
