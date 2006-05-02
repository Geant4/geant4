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
// $Id: G4TrajectoryModelFactories.cc,v 1.4 2006-05-02 20:47:40 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
