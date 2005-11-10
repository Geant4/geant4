//
// File name:     RadmonPhysicsPhotonEPDL.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsPhotonEPDL.cc,v 1.1 2005-11-10 08:15:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

#include "RadmonPhysicsPhotonEPDL.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"

#include "G4ProcessManager.hh"
#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyRayleigh.hh"
#include "G4StepLimiter.hh"

RadmonVSubPhysicsListWithLabel *                RadmonPhysicsPhotonEPDL :: New(void) const
{
 return new RadmonPhysicsPhotonEPDL;
}



void                                            RadmonPhysicsPhotonEPDL :: ConstructParticle(void)
{
 G4Gamma::GammaDefinition();
 G4Electron::ElectronDefinition();
}



void                                            RadmonPhysicsPhotonEPDL :: ConstructProcess(void)
{
 G4ProcessManager * manager(G4Gamma::GammaDefinition()->GetProcessManager());

 manager->AddDiscreteProcess(new G4LowEnergyPhotoElectric);
 manager->AddDiscreteProcess(new G4LowEnergyCompton);
 manager->AddDiscreteProcess(new G4LowEnergyGammaConversion);
 manager->AddDiscreteProcess(new G4LowEnergyRayleigh);
 manager->AddProcess(new G4StepLimiter(), ordInActive, ordInActive, 3);
}



void                                            RadmonPhysicsPhotonEPDL :: SetCuts(void)
{
}





const RadmonPhysicsInfoList &                   RadmonPhysicsPhotonEPDL :: Provides(void) const
{
 if (infoList.GetNPhysicsInfos()==0)
 {
  RadmonPhysicsInfo info;
  
  info.SetProcessName("PhotoElectric");
  info.SetParticleDefinition(G4Gamma::GammaDefinition());
  info.SetMinEnergy(250.*eV);
  info.SetMaxEnergy(100.*GeV);
  infoList.InsertPhysicsInfo(info);

  info.SetProcessName("Compton");
  infoList.InsertPhysicsInfo(info);

  info.SetProcessName("Rayleigh");
  infoList.InsertPhysicsInfo(info);

  info.SetProcessName("GammaConversion");
  info.SetMinEnergy(1.022000*MeV);
  info.SetMaxEnergy(100.*GeV);
  infoList.InsertPhysicsInfo(info);

  info.SetProcessName("StepLimiter");
  info.SetMinEnergy(0.*eV);
  info.SetMaxEnergy(DBL_MAX);
  infoList.InsertPhysicsInfo(info);
 }
 
 return infoList;
}
