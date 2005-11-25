//
// File name:     RadmonPhysicsDecay.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsDecay.cc,v 1.3 2005-11-25 11:53:02 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

#include "RadmonPhysicsDecay.hh"

#include "G4ProcessManager.hh"

#include "G4Decay.hh" 
#include "G4ParticleTable.hh"

RadmonVSubPhysicsListWithLabel *                RadmonPhysicsDecay :: New(void) const
{
 return new RadmonPhysicsDecay;
}



void                                            RadmonPhysicsDecay :: ConstructParticle(void)
{
}



void                                            RadmonPhysicsDecay :: ConstructProcess(void)
{
 G4ParticleTable * particleTable(G4ParticleTable::GetParticleTable());
 G4ParticleTable::G4PTblDicIterator & particleIterator(*particleTable->GetIterator());
    
 G4Decay * decayProcess(new G4Decay);

 particleIterator.reset();
 while (particleIterator())
 {
  G4ParticleDefinition * particle(particleIterator.value());
  G4ProcessManager * manager(particle->GetProcessManager());

  if (decayProcess->IsApplicable(*particle))
   manager->AddProcess(decayProcess, ordDefault, ordInActive, ordDefault);
 }
}



void                                            RadmonPhysicsDecay :: SetCuts(void)
{
}





const RadmonPhysicsInfoList &                   RadmonPhysicsDecay :: Provides(void) const
{
 if (infoList.GetNPhysicsInfos()==0)
 {
  G4ParticleTable * particleTable(G4ParticleTable::GetParticleTable());
  G4ParticleTable::G4PTblDicIterator & particleIterator(*particleTable->GetIterator());

  RadmonPhysicsInfo info;
  
  info.SetProcessName("AtRestDecay");
  info.SetMinEnergy(0*eV);
  info.SetMaxEnergy(1.e6*TeV);

  while (particleIterator())
  {  
   info.SetParticleDefinition(particleIterator.value());

   infoList.InsertPhysicsInfo(info);
  }
 }
 
 return infoList;
}
