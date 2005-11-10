//
// File name:     RadmonPhysicsNuclear.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsNuclear.cc,v 1.1 2005-11-10 08:15:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

#include "RadmonPhysicsNuclear.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

#include "G4ProcessManager.hh"

#include "G4ElectronNuclearProcess.hh"
#include "G4PositronNuclearProcess.hh"

#include "G4PhotoNuclearProcess.hh"

// models
#include "G4ElectroNuclearReaction.hh"

#include "G4GammaNuclearReaction.hh"
#include "G4TheoFSGenerator.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4GammaParticipants.hh"
#include "G4QGSModel.hh"


RadmonVSubPhysicsListWithLabel *                RadmonPhysicsNuclear :: New(void) const
{
 return new RadmonPhysicsNuclear;
}



void                                            RadmonPhysicsNuclear :: ConstructParticle(void)
{
 G4Positron::PositronDefinition();
 G4Gamma::GammaDefinition();
 G4Electron::ElectronDefinition();
}



void                                            RadmonPhysicsNuclear :: ConstructProcess(void)
{
 G4ElectroNuclearReaction* electronReaction(new G4ElectroNuclearReaction);
                                      
 G4ProcessManager * manager(G4Electron::ElectronDefinition()->GetProcessManager());
 G4ElectronNuclearProcess * electronNuclearProcess(new G4ElectronNuclearProcess);
 electronNuclearProcess->RegisterMe(electronReaction);
 manager->AddDiscreteProcess(electronNuclearProcess);

 manager=G4Positron::PositronDefinition()->GetProcessManager();
 G4PositronNuclearProcess * positronNuclearProcess(new G4PositronNuclearProcess);
 positronNuclearProcess->RegisterMe(electronReaction);
 manager->AddDiscreteProcess(positronNuclearProcess);


 G4GammaNuclearReaction * lowEGammaModel(new G4GammaNuclearReaction);
 lowEGammaModel->SetMaxEnergy(3.5*GeV);
     
 G4TheoFSGenerator * highEGammaModel(new G4TheoFSGenerator);
 highEGammaModel->SetTransport(new G4GeneratorPrecompoundInterface);
                                       
 G4QGSModel<G4GammaParticipants> * stringModel(new G4QGSModel<G4GammaParticipants>);
 stringModel->SetFragmentationModel(new G4ExcitedStringDecay(new G4QGSMFragmentation));
 
 highEGammaModel->SetHighEnergyGenerator(stringModel);
 highEGammaModel->SetMinEnergy(3.*GeV);
 highEGammaModel->SetMaxEnergy(100.*TeV); 
 
 manager=G4Gamma::GammaDefinition()->GetProcessManager();
 G4PhotoNuclearProcess * photoNuclearProcess(new G4PhotoNuclearProcess);
 photoNuclearProcess->RegisterMe(lowEGammaModel);
 photoNuclearProcess->RegisterMe(highEGammaModel);
 manager->AddDiscreteProcess(photoNuclearProcess);
}



void                                            RadmonPhysicsNuclear :: SetCuts(void)
{
}





const RadmonPhysicsInfoList &                   RadmonPhysicsNuclear :: Provides(void) const
{
 if (infoList.GetNPhysicsInfos()==0)
 {
  RadmonPhysicsInfo info;
  
  info.SetProcessName("NuclearReaction");
  info.SetParticleDefinition(G4Electron::ElectronDefinition());
  info.SetMinEnergy(0.*eV);
  info.SetMaxEnergy(30.*TeV);
  infoList.InsertPhysicsInfo(info);

  info.SetParticleDefinition(G4Positron::PositronDefinition());
  infoList.InsertPhysicsInfo(info);

  info.SetParticleDefinition(G4Gamma::GammaDefinition());
  info.SetMaxEnergy(100.*TeV);
  infoList.InsertPhysicsInfo(info);
 }
 
 return infoList;
}
