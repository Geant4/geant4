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
// File name:     RadmonPhysicsDecay.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsDecay.cc,v 1.5 2006-06-28 13:56:05 gunter Exp $
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

  particleIterator.reset();
  while (particleIterator())
  {  
   info.SetParticleDefinition(particleIterator.value());

   infoList.InsertPhysicsInfo(info);
  }
 }
 
 return infoList;
}
