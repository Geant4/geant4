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
// $Id: G4TrajectoryModelFactories.cc,v 1.2 2005/11/23 05:19:23 tinslay Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Jane Tinslay, John Allison, Joseph Perl October 2005

#include "G4ModelCommandsDrawByCharge.hh"
#include "G4ModelCommandsDrawByParticleID.hh"
#include "G4TrajectoryDrawByCharge.hh"
#include "G4TrajectoryDrawByParticleID.hh"
#include "G4TrajectoryModelFactories.hh"

G4TrajectoryDrawByChargeFactory::G4TrajectoryDrawByChargeFactory()
  :G4VModelFactory<G4VTrajectoryModel>("drawByCharge") 
{}

G4TrajectoryDrawByChargeFactory::~G4TrajectoryDrawByChargeFactory() {}

ModelAndMessengers
G4TrajectoryDrawByChargeFactory::Create(const G4String& placement, const G4String& name)
{
  // Create model
  G4TrajectoryDrawByCharge* model = new G4TrajectoryDrawByCharge(name);
  
  // Create associated messengers
  Messengers messengers;

  messengers.push_back(new G4ModelCommandDrawByChargeSet(model, placement));
  messengers.push_back(new G4ModelCommandDrawByChargeSetRGBA(model, placement));
  
  return ModelAndMessengers(model, messengers);
}

//Draw by particle ID
G4TrajectoryDrawByParticleIDFactory::G4TrajectoryDrawByParticleIDFactory()
  :G4VModelFactory<G4VTrajectoryModel>("drawByParticleID") 
{}

G4TrajectoryDrawByParticleIDFactory::~G4TrajectoryDrawByParticleIDFactory() {}

ModelAndMessengers
G4TrajectoryDrawByParticleIDFactory::Create(const G4String& placement, const G4String& name)
{
  // Create model
  G4TrajectoryDrawByParticleID* model = new G4TrajectoryDrawByParticleID(name);
  
  // Create associated messengers
  Messengers messengers;

  messengers.push_back(new G4ModelCommandDrawByParticleIDSet(model, placement));
  messengers.push_back(new G4ModelCommandDrawByParticleIDSetRGBA(model, placement));
  messengers.push_back(new G4ModelCommandDrawByParticleIDSetDefault(model, placement));
  messengers.push_back(new G4ModelCommandDrawByParticleIDSetDefaultRGBA(model, placement));
  
  return ModelAndMessengers(model, messengers);
}
