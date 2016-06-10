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
// $Id: G4TrajectoryParticleFilter.cc 66373 2012-12-18 09:41:34Z gcosmo $
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
G4TrajectoryParticleFilter::Evaluate(const G4VTrajectory& traj) const
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
