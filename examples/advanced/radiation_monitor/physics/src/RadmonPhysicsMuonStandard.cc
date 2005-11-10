//
// File name:     RadmonPhysicsMuonStandard.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsMuonStandard.cc,v 1.1 2005-11-10 08:15:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

#include "RadmonPhysicsMuonStandard.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"

#include "G4ProcessManager.hh"

#include "G4MultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuonMinusCaptureAtRest.hh" 
#include "G4StepLimiter.hh" 

RadmonVSubPhysicsListWithLabel *                RadmonPhysicsMuonStandard :: New(void) const
{
 return new RadmonPhysicsMuonStandard;
}



void                                            RadmonPhysicsMuonStandard :: ConstructParticle(void)
{
 G4Gamma::GammaDefinition();
 G4MuonPlus::MuonPlusDefinition();
 G4MuonMinus::MuonMinusDefinition();
 G4Electron::ElectronDefinition();
 G4Positron::PositronDefinition();
}



void                                            RadmonPhysicsMuonStandard :: ConstructProcess(void)
{
 G4ProcessManager * manager(G4MuonPlus::MuonPlusDefinition()->GetProcessManager());
 manager->AddProcess(new G4MultipleScattering, ordInActive,           1,           1);
 manager->AddProcess(new G4MuIonisation,       ordInActive,           2,           2);
 manager->AddProcess(new G4MuBremsstrahlung,   ordInActive,           3,           3);
 manager->AddProcess(new G4MuPairProduction,   ordInActive,           4,           4);
 manager->AddProcess(new G4StepLimiter,        ordInActive, ordInActive,           5);

 manager=G4MuonMinus::MuonMinusDefinition()->GetProcessManager();
 manager->AddProcess(new G4MultipleScattering, ordInActive,           1,           1);
 manager->AddProcess(new G4MuIonisation,       ordInActive,           2,           2);
 manager->AddProcess(new G4MuBremsstrahlung,   ordInActive,           3,           3);
 manager->AddProcess(new G4MuPairProduction,   ordInActive,           4,           4);
 manager->AddProcess(new G4MuonMinusCaptureAtRest,       0, ordInActive, ordInActive);
 manager->AddProcess(new G4StepLimiter,        ordInActive, ordInActive,           5);
}



void                                            RadmonPhysicsMuonStandard :: SetCuts(void)
{
}





const RadmonPhysicsInfoList &                   RadmonPhysicsMuonStandard :: Provides(void) const
{
 if (infoList.GetNPhysicsInfos()==0)
 {
  RadmonPhysicsInfo info;
  
  G4int i(2);
  
  info.SetParticleDefinition(G4MuonPlus::MuonPlusDefinition());
  while (i>0)
  {
   i--;
   info.SetProcessName("MultipleScattering");
   info.SetMinEnergy(0.1*keV);
   info.SetMaxEnergy(100.*GeV);
   infoList.InsertPhysicsInfo(info);

   info.SetProcessName("Ionisation");
   infoList.InsertPhysicsInfo(info);

   info.SetProcessName("Bremsstrahlung");
   infoList.InsertPhysicsInfo(info);

   info.SetProcessName("MuPairProduction");
   infoList.InsertPhysicsInfo(info);

   info.SetProcessName("StepLimiter");
   info.SetMinEnergy(0.*eV);
   info.SetMaxEnergy(DBL_MAX);
   infoList.InsertPhysicsInfo(info);
   
   info.SetParticleDefinition(G4MuonMinus::MuonMinusDefinition());
  }
  
  info.SetProcessName("CaptureAtRest");
  info.SetMinEnergy(0.*eV);
  info.SetMaxEnergy(0.*eV);
  infoList.InsertPhysicsInfo(info);
 }
 
 return infoList;
}
