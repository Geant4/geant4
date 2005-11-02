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
// $Id: G4TrajectoryDrawByParticleID.cc,v 1.1 2005-11-02 00:41:13 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay October 2005

#include "G4TrajectoryDrawByParticleID.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4TrajectoryDrawerUtils.hh"
#include "G4VisAttributes.hh"
#include "G4VTrajectory.hh"
#include "G4VVisManager.hh"

G4TrajectoryDrawByParticleID::G4TrajectoryDrawByParticleID(const G4String& name)
  :G4VTrajectoryDrawer(name)
  ,fDefault(G4Color::Grey)
{}

G4TrajectoryDrawByParticleID::~G4TrajectoryDrawByParticleID() {}

void
G4TrajectoryDrawByParticleID::SetDefault(const G4Colour& colour)
{
  fDefault = colour;
}

void
G4TrajectoryDrawByParticleID::Set(const G4String& particle, const G4String& colour)
{
  map<G4String, G4Colour>::iterator iter = fMap.find(particle);
  
  if (iter == fMap.end()) {
    G4Colour myColour(fDefault);

    if (!G4Colour::GetColour(colour, myColour)) {
      std::ostringstream o;
      o << "Failed to retrieve colour with key "<<colour<<" Using default " <<myColour;
      G4Exception
	("G4TrajectoryDrawByParticleID::Set(const G4String& particle, const G4String& colour)",
	 "MissingColourKey", JustWarning, o.str().c_str());
    }

    fMap[particle] = myColour;
  }
  else {
    std::ostringstream o;
    o << "Particle "<<particle<<" already has colour assigned ";
    G4Exception
      ("G4TrajectoryDrawByParticleID::Set(const G4String& particle, const G4String& colour)",
       "ParticleColourExists", JustWarning, o.str().c_str());
  }
}
/*
void
G4TrajectoryDrawByParticleID::Set(const G4String& particle, const G4Colour& colour)
{
  map<G4String, G4Colour>::iterator iter = fMap.find(particle);
  
  if (iter == fMap.end()) fMap[particle] = colour;
  else {
    std::ostringstream o;
    o << "Particle "<<particle<<" already has colour assigned ";
    G4Exception
      ("G4TrajectoryDrawByParticleID::Set(const G4String& particle, const G4Colour& colour)",
       "ParticleColourExists", JustWarning, o.str().c_str());
  }
}
*/
void
G4TrajectoryDrawByParticleID::Draw(const G4VTrajectory& traj, G4int i_mode) 
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (!pVVisManager) return;

  const G4double markerSize = std::abs(i_mode)/1000;
  G4bool lineRequired (i_mode >= 0);
  G4bool markersRequired (markerSize > 0.);   

  //Return if don't need to do anything
  if (!lineRequired && !markersRequired) return;

  //Get points to draw
  G4Polyline trajectoryLine;
  G4Polymarker stepPoints;
  G4Polymarker auxiliaryPoints;
  
  G4TrajectoryDrawerUtils::GetPoints(traj, trajectoryLine, 
				     auxiliaryPoints, stepPoints);
  
  if (lineRequired) {
    G4Colour colour(fDefault);
    G4String particle = traj.GetParticleName();

    map<G4String, G4Colour>::iterator iter = fMap.find(particle);
    if (iter != fMap.end()) colour = fMap[particle];

    G4VisAttributes trajectoryLineAttribs(colour);
    trajectoryLine.SetVisAttributes(&trajectoryLineAttribs);
    pVVisManager->Draw(trajectoryLine);
  }

  if (markersRequired) {
    auxiliaryPoints.SetMarkerType(G4Polymarker::squares);
    auxiliaryPoints.SetScreenSize(markerSize);
    auxiliaryPoints.SetFillStyle(G4VMarker::filled);
    G4VisAttributes auxiliaryPointsAttribs(G4Colour(0.,1.,1.));  // Magenta
    auxiliaryPoints.SetVisAttributes(&auxiliaryPointsAttribs);
    pVVisManager->Draw(auxiliaryPoints);
    
    stepPoints.SetMarkerType(G4Polymarker::circles);
    stepPoints.SetScreenSize(markerSize);
    stepPoints.SetFillStyle(G4VMarker::filled);  
    G4VisAttributes stepPointsAttribs(G4Colour(1.,1.,0.));  // Yellow.
    stepPoints.SetVisAttributes(&stepPointsAttribs);
    pVVisManager->Draw(stepPoints);
  }
  
  return;
}

void
G4TrajectoryDrawByParticleID::Print() const
{
  G4cout<<"G4TrajectoryDrawByParticleID drawer with name " <<GetName()<<G4endl;
  G4cout<<"Trajectory colour scheme: "<<G4endl;

  map<G4String, G4Colour>::const_iterator iter = fMap.begin();
  while (iter != fMap.end()) {
    G4cout<<"Particle, colour: "<<iter->first<<" : "<<iter->second<<G4endl;
    iter++;
  }
}
