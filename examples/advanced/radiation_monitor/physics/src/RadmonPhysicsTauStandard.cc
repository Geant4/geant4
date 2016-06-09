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
//
// File name:     RadmonPhysicsTauStandard.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsTauStandard.cc,v 1.3 2006/06/29 16:19:56 gunter Exp $
// Tag:           $Name: geant4-09-01 $
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
