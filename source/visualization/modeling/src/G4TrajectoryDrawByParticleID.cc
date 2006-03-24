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
// $Id: G4TrajectoryDrawByParticleID.cc,v 1.6 2006-03-24 20:22:43 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay, John Allison, Joseph Perl November 2005

#include "G4TrajectoryDrawByParticleID.hh"
#include "G4TrajectoryDrawerUtils.hh"
#include "G4VTrajectory.hh"
#include <sstream>

G4TrajectoryDrawByParticleID::G4TrajectoryDrawByParticleID(const G4String& name)
  :G4VTrajectoryModel(name)
  ,fDefault(G4Colour::Grey())
{}

G4TrajectoryDrawByParticleID::~G4TrajectoryDrawByParticleID() {}

void
G4TrajectoryDrawByParticleID::Draw(const G4VTrajectory& traj, const G4int& i_mode, const G4bool& visible) const
{
  G4Colour colour(fDefault);
  G4String particle = traj.GetParticleName();

  fMap.GetColour(particle, colour);

  G4TrajectoryDrawerUtils::DrawLineAndPoints(traj, i_mode, colour, visible);
}

void
G4TrajectoryDrawByParticleID::SetDefault(const G4String& colour)
{
  G4Colour myColour(G4Colour::White());      

  // Will not modify myColour if colour key does not exist  
  if (!G4Colour::GetColour(colour, myColour)) {
    std::ostringstream o;
    o << "G4Colour with key "<<colour<<" does not exist ";
    G4Exception
      ("G4TrajectoryDrawByParticleID::SetDefault(const G4String& colour)",
       "NonExistentColour", JustWarning, o.str().c_str());
  }

  SetDefault(myColour);
}

void
G4TrajectoryDrawByParticleID::SetDefault(const G4Colour& colour)
{
  fDefault = colour;
}

void
G4TrajectoryDrawByParticleID::Set(const G4String& particle, const G4String& colour)
{
  fMap.Set(particle, colour);
}

void
G4TrajectoryDrawByParticleID::Set(const G4String& particle, const G4Colour& colour)
{
  fMap[particle] = colour;
}

void
G4TrajectoryDrawByParticleID::Print(std::ostream& ostr) const
{
  ostr<<"G4TrajectoryDrawByParticleID model "<< Name() <<" colour scheme: "<<std::endl;
  
  ostr<<"Default colour: "<<fDefault<<G4endl;

  fMap.Print(ostr);
}
