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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPFSFissionFS.hh"
#include "G4ParticleHPManager.hh"
#include "G4ReactionProduct.hh"
#include "G4Nucleus.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"
#include "G4ThreeVector.hh"
#include "G4Poisson.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleHPDataUsed.hh"

void G4ParticleHPFSFissionFS::Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String &, G4ParticleDefinition*  )
{
  G4String tString = "/FS/";
  G4bool dbool;
  G4ParticleHPDataUsed aFile = theNames.GetName(static_cast<G4int>(A), static_cast<G4int>(Z), M, dirName, tString, dbool);
  G4String filename = aFile.GetName();
  SetAZMs( A, Z, M, aFile ); 
  if(!dbool)
  {
    hasAnyData = false;
    hasFSData = false; 
    hasXsec = false;
    return;
  }

  std::istringstream theData(std::ios::in);
  G4ParticleHPManager::GetInstance()->GetDataStream(filename,theData);

  G4int infoType, dataType;
  hasFSData = false; 
  while (theData >> infoType) // Loop checking, 11.05.2015, T. Koi
  {
    hasFSData = true; 
    theData >> dataType;
    switch(infoType)
    {
      case 1: 
        if(dataType==4) theNeutronAngularDis.Init(theData); 
        if(dataType==5) thePromptNeutronEnDis.Init(theData); 
        if(dataType==12) theFinalStatePhotons.InitMean(theData); 
        if(dataType==14) theFinalStatePhotons.InitAngular(theData); 
        if(dataType==15) theFinalStatePhotons.InitEnergies(theData); 
        break;
      case 2:
        if(dataType==1) theFinalStateNeutrons.InitMean(theData); 
        break;
      case 3:
        if(dataType==1) theFinalStateNeutrons.InitDelayed(theData); 
        if(dataType==5) theDelayedNeutronEnDis.Init(theData);
        break;
      case 4:
        if(dataType==1) theFinalStateNeutrons.InitPrompt(theData); 
        break;
      case 5:
        if(dataType==1) theEnergyRelease.Init(theData); 
        break;
      default:
        G4cout << "G4ParticleHPFSFissionFS::Init: unknown data type"<<dataType<<G4endl;
        throw G4HadronicException(__FILE__, __LINE__, "G4ParticleHPFSFissionFS::Init: unknown data type");
        break;
    }
  }
}

G4DynamicParticleVector*
G4ParticleHPFSFissionFS::ApplyYourself(G4int nPrompt, G4int nDelayed,
                                        G4double * theDecayConst)
{  
  G4int i;
  G4DynamicParticleVector * aResult = new G4DynamicParticleVector;
  G4ReactionProduct boosted;
  boosted.Lorentz( *(fCache.Get().theNeutronRP) , *(fCache.Get().theTarget) );
  G4double eKinetic = boosted.GetKineticEnergy();
    
  // Build neutrons
  G4ReactionProduct * theNeutrons = new G4ReactionProduct[nPrompt+nDelayed];
  for(i=0; i<nPrompt+nDelayed; ++i)
  {
    theNeutrons[i].SetDefinition(G4Neutron::Neutron());
  }
    
  // sample energies
  G4int it, dummy;
  G4double tempE;
  for(i=0; i<nPrompt; ++i)
  {
    tempE = thePromptNeutronEnDis.Sample(eKinetic, dummy); // energy distribution (file5) always in lab
    theNeutrons[i].SetKineticEnergy(tempE);
  }
  for(i=nPrompt; i<nPrompt+nDelayed; ++i)
  {
    theNeutrons[i].SetKineticEnergy(theDelayedNeutronEnDis.Sample(eKinetic, it));  // dito
    if(it==0) theNeutrons[i].SetKineticEnergy(thePromptNeutronEnDis.Sample(eKinetic, dummy));
    theDecayConst[i-nPrompt] = theFinalStateNeutrons.GetDecayConstant(it); // this is returned
  }

  // sample neutron angular distribution
  for(i=0; i<nPrompt+nDelayed; ++i)
  {
    theNeutronAngularDis.SampleAndUpdate(theNeutrons[i]); // angular comes back in lab automatically
  }
   
  // already in lab. Add neutrons to dynamic particle vector
  for(i=0; i<nPrompt+nDelayed; ++i)
  {
    G4DynamicParticle * dp = new G4DynamicParticle;
    dp->SetDefinition(theNeutrons[i].GetDefinition());
    dp->SetMomentum(theNeutrons[i].GetMomentum());
    aResult->push_back(dp);
  }
  delete [] theNeutrons;

  return aResult;
}

void G4ParticleHPFSFissionFS::SampleNeutronMult(G4int&all, G4int&Prompt, G4int&delayed, G4double eKinetic, G4int off)
{
   G4double promptNeutronMulti = 0;
   promptNeutronMulti = theFinalStateNeutrons.GetPrompt(eKinetic);
   G4double delayedNeutronMulti = 0;
   delayedNeutronMulti = theFinalStateNeutrons.GetDelayed(eKinetic);
   
   if(delayedNeutronMulti==0&&promptNeutronMulti==0)
   {
     Prompt = 0;
     delayed = 0;
     G4double totalNeutronMulti = theFinalStateNeutrons.GetMean(eKinetic);
     all = (G4int)G4Poisson(totalNeutronMulti-off);
     all += off;
   }
   else
   {   
     Prompt  = (G4int)G4Poisson(promptNeutronMulti-off);
     Prompt += off;
     delayed = (G4int)G4Poisson(delayedNeutronMulti);
     all = Prompt+delayed;
   }
}

G4DynamicParticleVector * G4ParticleHPFSFissionFS::GetPhotons()
{
   // sample photons
   G4ReactionProductVector * temp;
   G4ReactionProduct boosted;

   // the photon distributions are in the Nucleus rest frame.
   boosted.Lorentz( *(fCache.Get().theNeutronRP) , *(fCache.Get().theTarget) );
   G4double anEnergy = boosted.GetKineticEnergy();
   temp = theFinalStatePhotons.GetPhotons(anEnergy);
   if(temp == 0) { return 0; }

   // lorentz transform, and add photons to final state
   unsigned int i;
   G4DynamicParticleVector * result = new G4DynamicParticleVector;
   for(i=0; i<temp->size(); ++i)
   {
     // back to lab
     temp->operator[](i)->Lorentz(*(temp->operator[](i)), -1.* (*(fCache.Get().theTarget)) );
     G4DynamicParticle * theOne = new G4DynamicParticle;
     theOne->SetDefinition(temp->operator[](i)->GetDefinition());
     theOne->SetMomentum(temp->operator[](i)->GetMomentum());
     result->push_back(theOne);
     delete temp->operator[](i);
   }
   delete temp;
   return result;
}
