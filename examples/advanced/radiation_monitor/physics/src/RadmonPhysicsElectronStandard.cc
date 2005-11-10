//
// File name:     RadmonPhysicsElectronStandard.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsElectronStandard.cc,v 1.1 2005-11-10 08:15:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

#include "RadmonPhysicsElectronStandard.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"

#include "G4ProcessManager.hh"
#include "G4MultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4StepLimiter.hh"

RadmonVSubPhysicsListWithLabel *                RadmonPhysicsElectronStandard :: New(void) const
{
 return new RadmonPhysicsElectronStandard;
}



void                                            RadmonPhysicsElectronStandard :: ConstructParticle(void)
{
 G4Gamma::GammaDefinition();
 G4Electron::ElectronDefinition();
}



void                                            RadmonPhysicsElectronStandard :: ConstructProcess(void)
{
 G4ProcessManager * manager(G4Electron::ElectronDefinition()->GetProcessManager());

 manager->AddProcess(new G4MultipleScattering, ordInActive,           1, 1);
 manager->AddProcess(new G4eIonisation,        ordInActive,           2, 2);
 manager->AddProcess(new G4eBremsstrahlung,    ordInActive, ordInActive, 3);
 manager->AddProcess(new G4StepLimiter,        ordInActive, ordInActive, 4);
}



void                                            RadmonPhysicsElectronStandard :: SetCuts(void)
{
}





const RadmonPhysicsInfoList &                   RadmonPhysicsElectronStandard :: Provides(void) const
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
