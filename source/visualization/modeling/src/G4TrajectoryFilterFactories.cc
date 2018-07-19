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
/// $Id: G4TrajectoryFilterFactories.cc 98766 2016-08-09 14:17:17Z gcosmo $
//
//
// Trajectory filter model factories creating filters
// and associated messengers.
//
// Jane Tinslay March 2006
//
#include "G4AttributeFilterT.hh"
#include "G4ModelCommandsT.hh"
#include "G4TrajectoryFilterFactories.hh"
#include "G4TrajectoryChargeFilter.hh"
#include "G4TrajectoryParticleFilter.hh"
#include "G4TrajectoryOriginVolumeFilter.hh"
#include "G4TrajectoryEncounteredVolumeFilter.hh"

// Attribute filter
G4TrajectoryAttributeFilterFactory::G4TrajectoryAttributeFilterFactory()
  :G4VModelFactory< G4VFilter<G4VTrajectory> >("attributeFilter") 
{}

G4TrajectoryAttributeFilterFactory::~G4TrajectoryAttributeFilterFactory() {}

G4TrajectoryAttributeFilterFactory::ModelAndMessengers
G4TrajectoryAttributeFilterFactory::Create(const G4String& placement, const G4String& name)
{
  typedef G4AttributeFilterT<G4VTrajectory> G4TrajectoryAttributeFilter;
  // Create model
  G4TrajectoryAttributeFilter* model = new G4TrajectoryAttributeFilter(name);
  
  // Create associated messengers
  Messengers messengers;
  
  messengers.push_back(new G4ModelCmdSetString<G4TrajectoryAttributeFilter>(model, placement, "setAttribute"));
  messengers.push_back(new G4ModelCmdInvert<G4TrajectoryAttributeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdActive<G4TrajectoryAttributeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdVerbose<G4TrajectoryAttributeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdReset<G4TrajectoryAttributeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdAddInterval<G4TrajectoryAttributeFilter>(model, placement, "addInterval"));
  messengers.push_back(new G4ModelCmdAddValue<G4TrajectoryAttributeFilter>(model, placement, "addValue"));
 
  return ModelAndMessengers(model, messengers);
}

// Charge filter
G4TrajectoryChargeFilterFactory::G4TrajectoryChargeFilterFactory()
  :G4VModelFactory< G4VFilter<G4VTrajectory> >("chargeFilter") 
{}

G4TrajectoryChargeFilterFactory::~G4TrajectoryChargeFilterFactory() {}

G4TrajectoryChargeFilterFactory::ModelAndMessengers
G4TrajectoryChargeFilterFactory::Create(const G4String& placement, const G4String& name)
{
  // Create model
  G4TrajectoryChargeFilter* model = new G4TrajectoryChargeFilter(name);
  
  // Create associated messengers
  Messengers messengers;
  
  messengers.push_back(new G4ModelCmdAddString<G4TrajectoryChargeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdInvert<G4TrajectoryChargeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdActive<G4TrajectoryChargeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdVerbose<G4TrajectoryChargeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdReset<G4TrajectoryChargeFilter>(model, placement));
  
  return ModelAndMessengers(model, messengers);
}

// Particle type filter
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

// Origin volume filter
G4TrajectoryOriginVolumeFilterFactory::G4TrajectoryOriginVolumeFilterFactory()
:G4VModelFactory< G4VFilter<G4VTrajectory> >("originVolumeFilter")
{}

G4TrajectoryOriginVolumeFilterFactory::~G4TrajectoryOriginVolumeFilterFactory() {}

G4TrajectoryOriginVolumeFilterFactory::ModelAndMessengers
G4TrajectoryOriginVolumeFilterFactory::Create(const G4String& placement, const G4String& name)
{
  // Create model
  G4TrajectoryOriginVolumeFilter* model = new G4TrajectoryOriginVolumeFilter(name);

  // Create associated messengers
  Messengers messengers;

  messengers.push_back(new G4ModelCmdAddString<G4TrajectoryOriginVolumeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdInvert<G4TrajectoryOriginVolumeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdActive<G4TrajectoryOriginVolumeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdVerbose<G4TrajectoryOriginVolumeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdReset<G4TrajectoryOriginVolumeFilter>(model, placement));

  return ModelAndMessengers(model, messengers);
}

// Encountered volume filter
G4TrajectoryEncounteredVolumeFilterFactory::G4TrajectoryEncounteredVolumeFilterFactory()
:G4VModelFactory< G4VFilter<G4VTrajectory> >("encounteredVolumeFilter")
{}

G4TrajectoryEncounteredVolumeFilterFactory::~G4TrajectoryEncounteredVolumeFilterFactory() {}

G4TrajectoryEncounteredVolumeFilterFactory::ModelAndMessengers
G4TrajectoryEncounteredVolumeFilterFactory::Create(const G4String& placement, const G4String& name)
{
  // Create model
  G4TrajectoryEncounteredVolumeFilter* model = new G4TrajectoryEncounteredVolumeFilter(name);

  // Create associated messengers
  Messengers messengers;

  messengers.push_back(new G4ModelCmdAddString<G4TrajectoryEncounteredVolumeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdInvert<G4TrajectoryEncounteredVolumeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdActive<G4TrajectoryEncounteredVolumeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdVerbose<G4TrajectoryEncounteredVolumeFilter>(model, placement));
  messengers.push_back(new G4ModelCmdReset<G4TrajectoryEncounteredVolumeFilter>(model, placement));

  return ModelAndMessengers(model, messengers);
}
