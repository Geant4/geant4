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
// $Id: G4TrajectoryModelFactories.cc 98766 2016-08-09 14:17:17Z gcosmo $
//
// Jane Tinslay, John Allison, Joseph Perl October 2005

#include "G4ModelCompoundCommandsT.hh"
#include "G4ModelCommandsT.hh"
#include "G4ModelCommandUtils.hh"
#include "G4TrajectoryDrawByAttribute.hh"
#include "G4TrajectoryDrawByCharge.hh"
#include "G4TrajectoryDrawByOriginVolume.hh"
#include "G4TrajectoryDrawByParticleID.hh"
#include "G4TrajectoryDrawByEncounteredVolume.hh"
#include "G4TrajectoryGenericDrawer.hh"
#include "G4TrajectoryModelFactories.hh"
#include "G4VisTrajContext.hh"

// Draw by attribute
G4TrajectoryDrawByAttributeFactory::G4TrajectoryDrawByAttributeFactory()
  :G4VModelFactory<G4VTrajectoryModel>("drawByAttribute") 
{}

G4TrajectoryDrawByAttributeFactory::~G4TrajectoryDrawByAttributeFactory() {}

ModelAndMessengers
G4TrajectoryDrawByAttributeFactory::Create(const G4String& placement, const G4String& name)
{
  Messengers messengers;
  
  // Create default context and model 
  G4VisTrajContext* context = new G4VisTrajContext("default");

  G4TrajectoryDrawByAttribute* model = new G4TrajectoryDrawByAttribute(name, context);

  // Create messengers for default context configuration
  G4ModelCommandUtils::AddContextMsgrs(context, messengers, placement+"/"+name);

  messengers.push_back(new G4ModelCmdVerbose<G4TrajectoryDrawByAttribute>(model, placement));
  messengers.push_back(new G4ModelCmdSetString<G4TrajectoryDrawByAttribute>(model, placement, "setAttribute"));
  messengers.push_back(new G4ModelCmdAddIntervalContext<G4TrajectoryDrawByAttribute>(model, placement, "addInterval"));
  messengers.push_back(new G4ModelCmdAddValueContext<G4TrajectoryDrawByAttribute>(model, placement, "addValue"));

  return ModelAndMessengers(model, messengers);
}

G4TrajectoryGenericDrawerFactory::G4TrajectoryGenericDrawerFactory()
  :G4VModelFactory<G4VTrajectoryModel>("generic") 
{}

G4TrajectoryGenericDrawerFactory::~G4TrajectoryGenericDrawerFactory() {}

ModelAndMessengers
G4TrajectoryGenericDrawerFactory::Create(const G4String& placement, const G4String& name)
{
  Messengers messengers;
  
  // Create default context and model 
  G4VisTrajContext* context = new G4VisTrajContext("default");
  G4TrajectoryGenericDrawer* model = new G4TrajectoryGenericDrawer(name, context);

  // Create messengers for default context configuration
  G4ModelCommandUtils::AddContextMsgrs(context, messengers, placement+"/"+name);
  
  // Verbose command
  messengers.push_back(new G4ModelCmdVerbose<G4TrajectoryGenericDrawer>(model, placement));

  return ModelAndMessengers(model, messengers);
}

// drawByCharge
G4TrajectoryDrawByChargeFactory::G4TrajectoryDrawByChargeFactory()
  :G4VModelFactory<G4VTrajectoryModel>("drawByCharge") 
{}

G4TrajectoryDrawByChargeFactory::~G4TrajectoryDrawByChargeFactory() {}

ModelAndMessengers
G4TrajectoryDrawByChargeFactory::Create(const G4String& placement, const G4String& name)
{
  Messengers messengers;
  
  // Create default context and model 
  G4VisTrajContext* context = new G4VisTrajContext("default");
  G4TrajectoryDrawByCharge* model = new G4TrajectoryDrawByCharge(name, context);

  // Create messengers for default context configuration
  G4ModelCommandUtils::AddContextMsgrs(context, messengers, placement+"/"+name);
  
  // Create messengers for drawer
  messengers.push_back(new G4ModelCmdSetStringColour<G4TrajectoryDrawByCharge>(model, placement));
  messengers.push_back(new G4ModelCmdVerbose<G4TrajectoryDrawByCharge>(model, placement));

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
  Messengers messengers;

  // Create default context and model
  G4VisTrajContext* context = new G4VisTrajContext("default");
  G4TrajectoryDrawByParticleID* model = new G4TrajectoryDrawByParticleID(name, context);

  // Create messengers for default context configuration
  G4ModelCommandUtils::AddContextMsgrs(context, messengers, placement+"/"+name);

  // Create messengers for drawer
  messengers.push_back(new G4ModelCmdSetStringColour<G4TrajectoryDrawByParticleID>(model, placement));
  messengers.push_back(new G4ModelCmdSetDefaultColour<G4TrajectoryDrawByParticleID>(model, placement));
  messengers.push_back(new G4ModelCmdVerbose<G4TrajectoryDrawByParticleID>(model, placement));

  return ModelAndMessengers(model, messengers);
}

//Draw by origin volume
G4TrajectoryDrawByOriginVolumeFactory::G4TrajectoryDrawByOriginVolumeFactory()
:G4VModelFactory<G4VTrajectoryModel>("drawByOriginVolume")
{}

G4TrajectoryDrawByOriginVolumeFactory::~G4TrajectoryDrawByOriginVolumeFactory() {}

ModelAndMessengers
G4TrajectoryDrawByOriginVolumeFactory::Create(const G4String& placement, const G4String& name)
{
  Messengers messengers;

  // Create default context and model
  G4VisTrajContext* context = new G4VisTrajContext("default");
  G4TrajectoryDrawByOriginVolume* model = new G4TrajectoryDrawByOriginVolume(name, context);

  // Create messengers for default context configuration
  G4ModelCommandUtils::AddContextMsgrs(context, messengers, placement+"/"+name);

  // Create messengers for drawer
  messengers.push_back(new G4ModelCmdSetStringColour<G4TrajectoryDrawByOriginVolume>(model, placement));
  messengers.push_back(new G4ModelCmdSetDefaultColour<G4TrajectoryDrawByOriginVolume>(model, placement));
  messengers.push_back(new G4ModelCmdVerbose<G4TrajectoryDrawByOriginVolume>(model, placement));

  return ModelAndMessengers(model, messengers);
}

//Draw by encountered volume
G4TrajectoryDrawByEncounteredVolumeFactory::G4TrajectoryDrawByEncounteredVolumeFactory()
:G4VModelFactory<G4VTrajectoryModel>("drawByEncounteredVolume")
{}

G4TrajectoryDrawByEncounteredVolumeFactory::~G4TrajectoryDrawByEncounteredVolumeFactory() {}

ModelAndMessengers
G4TrajectoryDrawByEncounteredVolumeFactory::Create(const G4String& placement, const G4String& name)
{
  Messengers messengers;

  // Create default context and model
  G4VisTrajContext* context = new G4VisTrajContext("default");
  G4TrajectoryDrawByEncounteredVolume* model = new G4TrajectoryDrawByEncounteredVolume(name, context);

  // Create messengers for default context configuration
  G4ModelCommandUtils::AddContextMsgrs(context, messengers, placement+"/"+name);

  // Create messengers for drawer
  messengers.push_back(new G4ModelCmdSetStringColour<G4TrajectoryDrawByEncounteredVolume>(model, placement));
  messengers.push_back(new G4ModelCmdSetDefaultColour<G4TrajectoryDrawByEncounteredVolume>(model, placement));
  messengers.push_back(new G4ModelCmdVerbose<G4TrajectoryDrawByEncounteredVolume>(model, placement));

  return ModelAndMessengers(model, messengers);
}
