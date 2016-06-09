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
// File name:     RadmonPhysicsDecay.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsDecay.cc,v 1.6 2006/06/29 16:18:39 gunter Exp $
// Tag:           $Name: geant4-09-00 $
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
