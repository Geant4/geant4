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
/// $Id: G4TrajectoryFilterFactories.hh,v 1.2 2006-05-30 18:44:36 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Trajectory filter model factories creating filters
// and associated messengers.
//
// Jane Tinslay March 2006
//
#ifndef G4TRAJECTORYFILTERFACTORIES_HH
#define G4TRAJECTORYFILTERFACTORIES_HH

#include "G4VFilter.hh"
#include "G4VModelFactory.hh"
#include "G4VTrajectory.hh"

// Charge filter
class G4TrajectoryChargeFilterFactory : public G4VModelFactory< G4VFilter<G4VTrajectory>  > {

public: // With description

  typedef std::vector<G4UImessenger*> Messengers;
  typedef std::pair< G4VFilter<G4VTrajectory> *, Messengers > ModelAndMessengers;

  G4TrajectoryChargeFilterFactory();

  virtual ~G4TrajectoryChargeFilterFactory();
  
  ModelAndMessengers Create(const G4String& placement, const G4String& name);
    
};

// Particle filter
class G4TrajectoryParticleFilterFactory : public G4VModelFactory< G4VFilter<G4VTrajectory>  > {

public: // With description

  typedef std::vector<G4UImessenger*> Messengers;
  typedef std::pair< G4VFilter<G4VTrajectory> *, Messengers > ModelAndMessengers;

  G4TrajectoryParticleFilterFactory();

  virtual ~G4TrajectoryParticleFilterFactory();
  
  ModelAndMessengers Create(const G4String& placement, const G4String& name);
    
};

// Origin volume filter
class G4TrajectoryOriginVolumeFilterFactory : public G4VModelFactory< G4VFilter<G4VTrajectory>  > {

public: // With description

  typedef std::vector<G4UImessenger*> Messengers;
  typedef std::pair< G4VFilter<G4VTrajectory> *, Messengers > ModelAndMessengers;

  G4TrajectoryOriginVolumeFilterFactory();

  virtual ~G4TrajectoryOriginVolumeFilterFactory();
  
  ModelAndMessengers Create(const G4String& placement, const G4String& name);
    
};


#endif

