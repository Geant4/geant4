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
// File name:     RadmonPhysicsPhotonStandard.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsPhotonStandard.cc,v 1.2 2006-06-28 13:56:37 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//

#include "RadmonPhysicsPhotonStandard.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"

#include "G4ProcessManager.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4StepLimiter.hh"

RadmonVSubPhysicsListWithLabel *                RadmonPhysicsPhotonStandard :: New(void) const
{
 return new RadmonPhysicsPhotonStandard;
}



void                                            RadmonPhysicsPhotonStandard :: ConstructParticle(void)
{
 G4Gamma::GammaDefinition();
 G4Electron::ElectronDefinition();
}



void                                            RadmonPhysicsPhotonStandard :: ConstructProcess(void)
{
 G4ProcessManager * manager(G4Gamma::GammaDefinition()->GetProcessManager());

 manager->AddDiscreteProcess(new G4PhotoElectricEffect);
 manager->AddDiscreteProcess(new G4ComptonScattering);
 manager->AddDiscreteProcess(new G4GammaConversion);
 manager->AddProcess(new G4StepLimiter(), ordInActive, ordInActive, 4);
}



void                                            RadmonPhysicsPhotonStandard :: SetCuts(void)
{
}





const RadmonPhysicsInfoList &                   RadmonPhysicsPhotonStandard :: Provides(void) const
{
 if (infoList.GetNPhysicsInfos()==0)
 {
  RadmonPhysicsInfo info;
  
  info.SetProcessName("PhotoElectric");
  info.SetParticleDefinition(G4Gamma::GammaDefinition());
  info.SetMinEnergy(0.1*keV);
  info.SetMaxEnergy(100.*GeV);
  infoList.InsertPhysicsInfo(info);

  info.SetProcessName("Compton");
  infoList.InsertPhysicsInfo(info);

  info.SetProcessName("GammaConversion");
  info.SetMinEnergy(2.0*electron_mass_c2);
  info.SetMaxEnergy(100.*GeV);
  infoList.InsertPhysicsInfo(info);

  info.SetProcessName("StepLimiter");
  info.SetMinEnergy(0.*eV);
  info.SetMaxEnergy(DBL_MAX);
  infoList.InsertPhysicsInfo(info);
 }
 
 return infoList;
}
