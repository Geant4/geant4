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
// $Id: G4TrajectoryModelFactories.cc,v 1.5 2006/06/29 21:33:18 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// Jane Tinslay, John Allison, Joseph Perl October 2005

#include "G4ModelCommandsT.hh"
#include "G4TrajectoryDrawByCharge.hh"
#include "G4TrajectoryDrawByOriginVolume.hh"
#include "G4TrajectoryDrawByParticleID.hh"
#include "G4TrajectoryGenericDrawer.hh"
#include "G4TrajectoryModelFactories.hh"
#include "G4VisTrajContext.hh"

namespace {
 
  void AddContextMsgrs(G4VisTrajContext* context, std::vector<G4UImessenger*>& messengers,
		       const G4String& placement)
  {
    messengers.push_back(new G4ModelCmdCreateContextDir<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetDrawLine<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetLineVisible<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetLineColour<G4VisTrajContext>(context, placement));
    
    messengers.push_back(new G4ModelCmdSetDrawStepPts<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetStepPtsVisible<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetStepPtsColour<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetStepPtsSize<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetStepPtsType<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetStepPtsFillStyle<G4VisTrajContext>(context, placement));

    messengers.push_back(new G4ModelCmdSetDrawAuxPts<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetAuxPtsVisible<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetAuxPtsColour<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetAuxPtsSize<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetAuxPtsType<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetAuxPtsFillStyle<G4VisTrajContext>(context, placement));
  }
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
  AddContextMsgrs(context, messengers, placement+"/"+name);
  
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
  AddContextMsgrs(context, messengers, placement+"/"+name);
  
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
  AddContextMsgrs(context, messengers, placement+"/"+name);

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
  AddContextMsgrs(context, messengers, placement+"/"+name);

  // Create messengers for drawer
  messengers.push_back(new G4ModelCmdSetStringColour<G4TrajectoryDrawByOriginVolume>(model, placement));
  messengers.push_back(new G4ModelCmdSetDefaultColour<G4TrajectoryDrawByOriginVolume>(model, placement));
  messengers.push_back(new G4ModelCmdVerbose<G4TrajectoryDrawByOriginVolume>(model, placement));

  return ModelAndMessengers(model, messengers);
}
