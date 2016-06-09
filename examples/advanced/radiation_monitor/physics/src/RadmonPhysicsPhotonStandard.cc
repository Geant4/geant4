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
// File name:     RadmonPhysicsPhotonStandard.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsPhotonStandard.cc,v 1.3 2006/06/29 16:19:48 gunter Exp $
// Tag:           $Name: geant4-09-01 $
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
