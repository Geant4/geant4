//
// File name:     RadmonPhysicsElectronEEDL.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsElectronEEDL.cc,v 1.2 2005-12-06 19:37:38 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

#include "RadmonPhysicsElectronEEDL.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4LeptonConstructor.hh"

#include "G4ProcessManager.hh"
#include "G4MultipleScattering.hh"
#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4StepLimiter.hh"

RadmonVSubPhysicsListWithLabel *                RadmonPhysicsElectronEEDL :: New(void) const
{
 return new RadmonPhysicsElectronEEDL;
}



void                                            RadmonPhysicsElectronEEDL :: ConstructParticle(void)
{
 G4Gamma::GammaDefinition();
 G4LeptonConstructor::ConstructParticle();
}



void                                            RadmonPhysicsElectronEEDL :: ConstructProcess(void)
{
 G4ProcessManager * manager(G4Electron::ElectronDefinition()->GetProcessManager());

 manager->AddProcess(new G4MultipleScattering,      ordInActive,           1, 1);
 manager->AddProcess(new G4LowEnergyIonisation,     ordInActive,           2, 2);
 manager->AddProcess(new G4LowEnergyBremsstrahlung, ordInActive, ordInActive, 3);
 manager->AddProcess(new G4StepLimiter,             ordInActive, ordInActive, 4);
}



void                                            RadmonPhysicsElectronEEDL :: SetCuts(void)
{
}





const RadmonPhysicsInfoList &                   RadmonPhysicsElectronEEDL :: Provides(void) const
{
 if (infoList.GetNPhysicsInfos()==0)
 {
  RadmonPhysicsInfo info;
  
  info.SetProcessName("MultipleScattering");
  info.SetParticleDefinition(G4Electron::ElectronDefinition());
  info.SetMinEnergy(0.1*keV);
  info.SetMaxEnergy(100.*TeV);
  infoList.InsertPhysicsInfo(info);

  info.SetProcessName("Ionisation");
  info.SetMinEnergy(250.*eV);
  info.SetMaxEnergy(100.*GeV);
  infoList.InsertPhysicsInfo(info);

  info.SetProcessName("Bremsstrahlung");
  infoList.InsertPhysicsInfo(info);

  info.SetProcessName("StepLimiter");
  info.SetMinEnergy(0.*eV);
  info.SetMaxEnergy(DBL_MAX);
  infoList.InsertPhysicsInfo(info);
 }
 
 return infoList;
}
