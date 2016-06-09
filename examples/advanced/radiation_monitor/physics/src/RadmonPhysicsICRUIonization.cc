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
// File name:     RadmonPhysicsICRUIonization.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsICRUIonization.cc,v 1.5 2006/06/29 16:18:49 gunter Exp $
// Tag:           $Name: geant4-08-02 $
//

#include "RadmonPhysicsICRUIonization.hh"

#include "G4ProcessManager.hh"

#include "G4MultipleScattering.hh" 
#include "G4hIonisation.hh" 
#include "G4StepLimiter.hh" 

#include "G4ParticleTable.hh"

RadmonVSubPhysicsListWithLabel *                RadmonPhysicsICRUIonization :: New(void) const
{
 return new RadmonPhysicsICRUIonization;
}



void                                            RadmonPhysicsICRUIonization :: ConstructParticle(void)
{
}



void                                            RadmonPhysicsICRUIonization :: ConstructProcess(void)
{
 G4ParticleTable * particleTable(G4ParticleTable::GetParticleTable());
 G4ParticleTable::G4PTblDicIterator & particleIterator(*particleTable->GetIterator());
 
 particleIterator.reset();
 while (particleIterator())
 {
  G4ParticleDefinition * particle(particleIterator.value());
  G4String particleName(particle->GetParticleName());

  if ((particle->GetPDGCharge()!=0) && (!particle -> IsShortLived()) && 
      particleName!="e+" && particleName!="mu+" && particleName!="tau+" &&
      particleName!="e-" && particleName!="mu-" && particleName!="tau-" && particleName!="chargedgeantino")
  {
   G4ProcessManager * manager(particle->GetProcessManager());
 
   manager->AddProcess(new G4MultipleScattering(), ordInActive, 1,           1);
   manager->AddProcess(new G4hIonisation(),        ordInActive, 2,           2);
   manager->AddProcess(new G4StepLimiter(),        ordInActive, ordInActive, 3);
  }
 }
}



void                                            RadmonPhysicsICRUIonization :: SetCuts(void)
{
}





const RadmonPhysicsInfoList &                   RadmonPhysicsICRUIonization :: Provides(void) const
{
 if (infoList.GetNPhysicsInfos()==0)
 {
  G4ParticleTable * particleTable(G4ParticleTable::GetParticleTable());
  G4ParticleTable::G4PTblDicIterator & particleIterator(*particleTable->GetIterator());

  RadmonPhysicsInfo info[3];
  
  info[0].SetProcessName("MultipleScattering");
  info[0].SetMinEnergy(0*eV);
  info[0].SetMaxEnergy(1.e6*TeV);

  info[1].SetProcessName("Ionisation");
  info[1].SetMinEnergy(100.*eV);
  info[1].SetMaxEnergy(100.*TeV);

  info[2].SetProcessName("StepLimiter");
  info[2].SetMinEnergy(0*eV);
  info[2].SetMaxEnergy(DBL_MAX);

  particleIterator.reset();
  while (particleIterator())
  {
   G4ParticleDefinition * particle(particleIterator.value());
   G4String particleName(particle->GetParticleName());

   if ((particle->GetPDGCharge()!=0) && (!particle -> IsShortLived()) && 
        particleName!="e+" && particleName!="mu+" && particleName!="tau+" &&
        particleName!="e-" && particleName!="mu-" && particleName!="tau-" && particleName!="chargedgeantino")
   {
    G4int i(3);
   
    while (i>0)
    {
     i--;
    
     info[i].SetParticleDefinition(particleIterator.value());
     infoList.InsertPhysicsInfo(info[i]);
    }
   }
  }
 }
 
 return infoList;
}
