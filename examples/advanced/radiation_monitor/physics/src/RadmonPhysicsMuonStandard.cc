//
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
//
// File name:     RadmonPhysicsMuonStandard.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsMuonStandard.cc,v 1.3 2006-06-28 13:56:27 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//

#include "RadmonPhysicsMuonStandard.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4LeptonConstructor.hh"

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
 G4LeptonConstructor::ConstructParticle();
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
