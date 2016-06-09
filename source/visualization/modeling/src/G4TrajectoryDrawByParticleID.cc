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
// $Id: G4TrajectoryDrawByParticleID.cc,v 1.4 2005/11/23 20:24:15 tinslay Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Jane Tinslay, John Allison, Joseph Perl November 2005

#include "G4TrajectoryDrawByParticleID.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4TrajectoryDrawerUtils.hh"
#include "G4VisAttributes.hh"
#include "G4VTrajectory.hh"
#include "G4VVisManager.hh"
#include <sstream>

G4TrajectoryDrawByParticleID::G4TrajectoryDrawByParticleID(const G4String& name)
  :G4VTrajectoryModel(name)
  ,fDefault(G4Color::Grey())
{}

G4TrajectoryDrawByParticleID::~G4TrajectoryDrawByParticleID() {}

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
  G4Colour myColour(fDefault);

  // Will not modify myColour if colour key does not exist  
  
  if (!G4Colour::GetColour(colour, myColour)) {
    std::ostringstream o;
    o << "G4Colour with key "<<colour<<" does not exist ";
    G4Exception
      ("G4TrajectoryDrawByParticleID::Set(const G4String& particle, const G4String& colour)",
       "NonExistentColour", JustWarning, o.str().c_str());
  }

  Set(particle, myColour);
}

void
G4TrajectoryDrawByParticleID::Set(const G4String& particle, const G4Colour& colour)
{
  std::map<G4String, G4Colour>::iterator iter = fMap.find(particle);
  
  if (iter == fMap.end()) {
    fMap[particle] = colour;
  }
  else {
    std::ostringstream o;
    o << "Particle "<<particle<<" already has colour "<<colour<<" assigned.";
    G4Exception
      ("G4TrajectoryDrawByParticleID::Set(const G4String& particle, const G4String& colour)",
       "ParticleColourExists", JustWarning, o.str().c_str());
  }
}

void
G4TrajectoryDrawByParticleID::Draw(const G4VTrajectory& traj, G4int i_mode) const
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (0 == pVVisManager) return;

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

    map<G4String, G4Colour>::const_iterator iter = fMap.find(particle);
    if (iter != fMap.end()) colour = iter->second;

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
G4TrajectoryDrawByParticleID::Print(std::ostream& ostr) const
{
  ostr<<"G4TrajectoryDrawByParticleID model "<< Name() <<" colour scheme: "<<std::endl;
  
  ostr<<"Default : "<<fDefault<<G4endl;

  std::map<G4String, G4Colour>::const_iterator iter = fMap.begin();
  
  while (iter != fMap.end()) {
    ostr<<iter->first <<" : "<< iter->second <<G4endl;
    iter++;
  }
}
