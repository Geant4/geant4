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
// $Id: G4TrajectoryDrawByParticleID.cc,v 1.8 2006-06-29 21:33:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay, John Allison, Joseph Perl November 2005

#include "G4TrajectoryDrawByParticleID.hh"
#include "G4TrajectoryDrawerUtils.hh"
#include "G4VisTrajContext.hh"
#include "G4VTrajectory.hh"
#include <sstream>

G4TrajectoryDrawByParticleID::G4TrajectoryDrawByParticleID(const G4String& name, G4VisTrajContext* context)
  :G4VTrajectoryModel(name, context)
  ,fDefault(G4Colour::Grey())
{}

G4TrajectoryDrawByParticleID::~G4TrajectoryDrawByParticleID() {}

void
G4TrajectoryDrawByParticleID::Draw(const G4VTrajectory& traj, const G4int& i_mode, const G4bool& visible) const
{
  G4Colour colour(fDefault);
  G4String particle = traj.GetParticleName();

  fMap.GetColour(particle, colour);

  G4VisTrajContext myContext(GetContext());
  
  myContext.SetLineColour(colour);
  myContext.SetVisible(visible);
  
  if (GetVerbose()) {
    G4cout<<"G4TrajectoryDrawByParticleID drawer named "<<Name();
    G4cout<<", drawing trajectory with particle type, "<<particle<<G4endl;
    G4cout<<", with configuration:"<<G4endl;
    myContext.Print(G4cout);
  }

  G4TrajectoryDrawerUtils::DrawLineAndPoints(traj, myContext, i_mode);
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

  ostr<<"Default configuration:"<<G4endl;
  GetContext().Print(G4cout);
}
