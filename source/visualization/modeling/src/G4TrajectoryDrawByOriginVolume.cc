// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4TrajectoryDrawByOriginVolume.cc,v 1.2 2006-03-24 20:22:43 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay March 2006

#include "G4TrajectoryDrawByOriginVolume.hh"
#include "G4TrajectoryDrawerUtils.hh"
#include "G4VTrajectory.hh"
#include "G4TransportationManager.hh"
#include "G4VTrajectoryPoint.hh"
#include <sstream>

G4TrajectoryDrawByOriginVolume::G4TrajectoryDrawByOriginVolume(const G4String& name)
  :G4VTrajectoryModel(name)
  ,fDefault(G4Colour::Grey())
{}

G4TrajectoryDrawByOriginVolume::~G4TrajectoryDrawByOriginVolume() {}

void
G4TrajectoryDrawByOriginVolume::Draw(const G4VTrajectory& traj, const G4int& i_mode, const G4bool& visible) const
{
  G4Colour colour(fDefault);
  
  G4Navigator* navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  
  G4VTrajectoryPoint* aTrajectoryPoint = traj.GetPoint(0);
  assert (0 != aTrajectoryPoint);
  
  G4VPhysicalVolume* volume = navigator->LocateGlobalPointAndSetup(aTrajectoryPoint->GetPosition());

  // Logical volumes form basis.
  G4LogicalVolume* logicalVolume = volume->GetLogicalVolume();
  assert (0 != logicalVolume);

  G4String logicalName = logicalVolume->GetName();
  fMap.GetColour(logicalName, colour);

  // Override with physical volume colouring if it exists
  G4String physicalName = volume->GetName();
  fMap.GetColour(physicalName, colour);
   
  G4TrajectoryDrawerUtils::DrawLineAndPoints(traj, i_mode, colour, visible);
}

void
G4TrajectoryDrawByOriginVolume::SetDefault(const G4String& colour)
{
  G4Colour myColour;      

  // Will not modify default colour if colour key does not exist  
  if (!G4Colour::GetColour(colour, myColour)) {
    std::ostringstream o;
    o << "G4Colour with key "<<colour<<" does not exist ";
    G4Exception
      ("G4TrajectoryDrawByOriginParticleID::SetDefault(const G4String& colour)",
       "NonExistentColour", JustWarning, o.str().c_str());
    return;
  }

  SetDefault(myColour);
}

void
G4TrajectoryDrawByOriginVolume::SetDefault(const G4Colour& colour)
{
  fDefault = colour;
}

void
G4TrajectoryDrawByOriginVolume::Set(const G4String& particle, const G4String& colour)
{
  fMap.Set(particle, colour);
}

void
G4TrajectoryDrawByOriginVolume::Set(const G4String& particle, const G4Colour& colour)
{
  fMap[particle] = colour;
}

void
G4TrajectoryDrawByOriginVolume::Print(std::ostream& ostr) const
{
  ostr<<"G4TrajectoryDrawByOriginVolume model "<< Name() <<" colour scheme: "<<std::endl;
  ostr<<"Default : "<<fDefault<<G4endl;

  fMap.Print(ostr);
}
