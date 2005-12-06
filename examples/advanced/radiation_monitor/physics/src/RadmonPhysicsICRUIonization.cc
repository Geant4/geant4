//
// File name:     RadmonPhysicsICRUIonization.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsICRUIonization.cc,v 1.3 2005-12-06 19:37:17 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
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
