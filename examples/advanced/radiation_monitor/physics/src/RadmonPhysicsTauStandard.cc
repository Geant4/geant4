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
// File name:     RadmonPhysicsTauStandard.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsTauStandard.cc,v 1.2 2006-06-28 13:56:45 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//

#include "RadmonPhysicsTauStandard.hh"

#include "G4Electron.hh"
#include "G4TauPlus.hh"
#include "G4TauMinus.hh"

#include "G4ProcessManager.hh"

#include "G4MultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4StepLimiter.hh"

RadmonVSubPhysicsListWithLabel *                RadmonPhysicsTauStandard :: New(void) const
{
 return new RadmonPhysicsTauStandard;
}



void                                            RadmonPhysicsTauStandard :: ConstructParticle(void)
{
 G4TauPlus::TauPlusDefinition();
 G4TauMinus::TauMinusDefinition();
 G4Electron::ElectronDefinition();
}



void                                            RadmonPhysicsTauStandard :: ConstructProcess(void)
{
 G4ProcessManager * manager(G4TauPlus::TauPlusDefinition()->GetProcessManager());
 manager->AddProcess(new G4MultipleScattering, ordInActive,           1,           1);
 manager->AddProcess(new G4hIonisation,        ordInActive,           2,           2);
 manager->AddProcess(new G4StepLimiter,        ordInActive, ordInActive,           3);

 manager=G4TauMinus::TauMinusDefinition()->GetProcessManager();
 manager->AddProcess(new G4MultipleScattering, ordInActive,           1,           1);
 manager->AddProcess(new G4hIonisation,        ordInActive,           2,           2);
 manager->AddProcess(new G4StepLimiter,        ordInActive, ordInActive,           3);
}



void                                            RadmonPhysicsTauStandard :: SetCuts(void)
{
}





const RadmonPhysicsInfoList &                   RadmonPhysicsTauStandard :: Provides(void) const
{
 if (infoList.GetNPhysicsInfos()==0)
 {
  RadmonPhysicsInfo info;
  
  G4int i(2);
  
  info.SetParticleDefinition(G4TauPlus::TauPlusDefinition());
  while (i>0)
  {
   i--;
   info.SetProcessName("MultipleScattering");
   info.SetMinEnergy(0.1*keV);
   info.SetMaxEnergy(100.*GeV);
   infoList.InsertPhysicsInfo(info);

   info.SetProcessName("Ionisation");
   infoList.InsertPhysicsInfo(info);

   info.SetProcessName("StepLimiter");
   info.SetMinEnergy(0.*eV);
   info.SetMaxEnergy(DBL_MAX);
   infoList.InsertPhysicsInfo(info);
   
   info.SetParticleDefinition(G4TauMinus::TauMinusDefinition());
  }
 }
 
 return infoList;
}
