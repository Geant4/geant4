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
// $Id: G4TrajectoryParticleFilter.cc,v 1.1 2006-03-28 18:01:18 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Filter trajectories according to particle type. Only registered 
// particle types will pass the filter.
//
// Jane Tinslay March 2006
//
#include "G4TrajectoryParticleFilter.hh"

G4TrajectoryParticleFilter::G4TrajectoryParticleFilter(const G4String& name)
  :G4SmartFilter<G4VTrajectory>(name)
{}

G4TrajectoryParticleFilter::~G4TrajectoryParticleFilter() {}

bool
G4TrajectoryParticleFilter::Evaluate(const G4VTrajectory& traj)
{
  G4String particle = traj.GetParticleName();

  if (GetVerbose()) G4cout<<"G4TrajectoryParticleFilter processing trajectory with particle type: "<<particle<<G4endl;

  std::vector<G4String>::const_iterator iter = std::find(fParticles.begin(), fParticles.end(), particle);

  // Fail if particle type not found in particle list
  if (iter == fParticles.end()) return false;

  return true;
}


void
G4TrajectoryParticleFilter::Add(const G4String& particle)
{
  fParticles.push_back(particle);
}

void
G4TrajectoryParticleFilter::Print(std::ostream& ostr) const
{
  ostr<<"Particle types registered: "<<G4endl;
  std::vector<G4String>::const_iterator iter = fParticles.begin();
  
  while (iter != fParticles.end()) {
    ostr<<*iter<<G4endl;    
    iter++;
  }
}

void 
G4TrajectoryParticleFilter::Clear()
{
  // Clear particle type vector
  fParticles.clear();
}
