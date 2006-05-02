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
/// $Id: G4TrajectoryFilterFactories.cc,v 1.2 2006-05-02 20:47:40 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Trajectory filter model factories creating filters
// and associated messengers.
//
// Jane Tinslay March 2006
//
#include "G4ModelCommandsT.hh"
#include "G4TrajectoryFilterFactories.hh"
#include "G4TrajectoryParticleFilter.hh"

G4TrajectoryParticleFilterFactory::G4TrajectoryParticleFilterFactory()
  :G4VModelFactory< G4VFilter<G4VTrajectory> >("particleFilter") 
{}

G4TrajectoryParticleFilterFactory::~G4TrajectoryParticleFilterFactory() {}

G4TrajectoryParticleFilterFactory::ModelAndMessengers
G4TrajectoryParticleFilterFactory::Create(const G4String& placement, const G4String& name)
{
  // Create model
  G4TrajectoryParticleFilter* model = new G4TrajectoryParticleFilter(name);
  
  // Create associated messengers
  Messengers messengers;
  
  messengers.push_back(new G4ModelCmdAddString<G4TrajectoryParticleFilter>(model, placement));
  messengers.push_back(new G4ModelCmdInvert<G4TrajectoryParticleFilter>(model, placement));
  messengers.push_back(new G4ModelCmdActive<G4TrajectoryParticleFilter>(model, placement));
  messengers.push_back(new G4ModelCmdVerbose<G4TrajectoryParticleFilter>(model, placement));
  messengers.push_back(new G4ModelCmdReset<G4TrajectoryParticleFilter>(model, placement));
  
  return ModelAndMessengers(model, messengers);
}

