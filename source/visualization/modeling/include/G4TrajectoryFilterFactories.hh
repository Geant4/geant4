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
/// $Id: G4TrajectoryFilterFactories.hh 98766 2016-08-09 14:17:17Z gcosmo $
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

// Attribute filter
class G4TrajectoryAttributeFilterFactory : public G4VModelFactory< G4VFilter<G4VTrajectory>  > {

public: // With description

  typedef std::vector<G4UImessenger*> Messengers;
  typedef std::pair< G4VFilter<G4VTrajectory> *, Messengers > ModelAndMessengers;

  G4TrajectoryAttributeFilterFactory();

  virtual ~G4TrajectoryAttributeFilterFactory();
  
  ModelAndMessengers Create(const G4String& placement, const G4String& name);
    
};

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

// Encountered volume filter
class G4TrajectoryEncounteredVolumeFilterFactory : public G4VModelFactory< G4VFilter<G4VTrajectory>  > {

public: // With description

  typedef std::vector<G4UImessenger*> Messengers;
  typedef std::pair< G4VFilter<G4VTrajectory> *, Messengers > ModelAndMessengers;

  G4TrajectoryEncounteredVolumeFilterFactory();

  virtual ~G4TrajectoryEncounteredVolumeFilterFactory();

  ModelAndMessengers Create(const G4String& placement, const G4String& name);

};


#endif

