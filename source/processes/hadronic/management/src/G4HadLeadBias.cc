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
// 20110906  M. Kelsey -- Use reference to G4HadSecondary instead of pointer

#include "G4HadLeadBias.hh"
#include "G4Gamma.hh"
#include "G4PionZero.hh"
#include "Randomize.hh"
#include "G4HadFinalState.hh"

  G4HadFinalState * G4HadLeadBias::Bias(G4HadFinalState * result)
  {
    // G4cerr << "bias enter"<<G4endl;
    G4int nMeson(0), nBaryon(0), npi0(0), ngamma(0), nLepton(0);
    G4int i(0);
    G4int maxE = -1;
    G4double emax = 0;
    if(result->GetStatusChange()==isAlive) 
    {
      emax = result->GetEnergyChange();
    }
    //G4cout << "max energy "<<G4endl;
    for(i=0;i<result->GetNumberOfSecondaries();i++)
    {
      if(result->GetSecondary(i)->GetParticle()->GetKineticEnergy()>emax)
      {
        maxE = i;
	emax = result->GetSecondary(i)->GetParticle()->GetKineticEnergy();
      }
    }
    //G4cout <<"loop1"<<G4endl;
    for(i=0; i<result->GetNumberOfSecondaries(); i++)
    {
      const G4DynamicParticle* aSecTrack = result->GetSecondary(i)->GetParticle();
      if(i==maxE)
      {
      }
      else if(aSecTrack->GetDefinition()->GetBaryonNumber()!=0) 
      {
        nBaryon++;
      }
      else if(aSecTrack->GetDefinition()->GetLeptonNumber()!=0) 
      {
        nLepton++;
      }
      else if(aSecTrack->GetDefinition()==G4Gamma::Gamma())
      {
        ngamma++;
      }
      else if(aSecTrack->GetDefinition()==G4PionZero::PionZero())
      {
        npi0++;
      }
      else
      {
        nMeson++;
      }
    }
     //G4cout << "BiasDebug 1 = "<<result->GetNumberOfSecondaries()<<" "
     //       <<nMeson<<" "<< nBaryon<<" "<< npi0<<" "<< ngamma<<" "<< nLepton<<G4endl;
    G4double mesonWeight = nMeson;
    G4double baryonWeight = nBaryon;
    G4double gammaWeight = ngamma;
    G4double npi0Weight = npi0;
    G4double leptonWeight = nLepton;
    G4int randomMeson = static_cast<G4int>((nMeson+1)*G4UniformRand());
    G4int randomBaryon = static_cast<G4int>((nBaryon+1)*G4UniformRand());
    G4int randomGamma = static_cast<G4int>((ngamma+1)*G4UniformRand());
    G4int randomPi0 = static_cast<G4int>((npi0+1)*G4UniformRand());
    G4int randomLepton = static_cast<G4int>((nLepton+1)*G4UniformRand());
    
    std::vector<G4HadSecondary> buffer;
    G4int cMeson(0), cBaryon(0), cpi0(0), cgamma(0), cLepton(0);
    for(i=0; i<result->GetNumberOfSecondaries(); i++)
    {
      G4bool aCatch = false;
      G4double weight = 1;
      const G4HadSecondary* aSecTrack = result->GetSecondary(i);
      G4ParticleDefinition* aSecDef = aSecTrack->GetParticle()->GetDefinition();
      if(i==maxE)
      {
        aCatch = true;
	weight = 1;
      }
      else if(aSecDef->GetBaryonNumber()!=0) 
      {
	if(++cBaryon==randomBaryon) 
	{
	  aCatch = true;
	  weight = baryonWeight;
	}
      }
      else if(aSecDef->GetLeptonNumber()!=0) 
      {
        if(++cLepton==randomLepton) 
	{
	  aCatch = true;
	  weight = leptonWeight;
	}
      }
      else if(aSecDef==G4Gamma::Gamma())
      {
        if(++cgamma==randomGamma) 
	{
	  aCatch = true;
	  weight = gammaWeight;
	}
      }
      else if(aSecDef==G4PionZero::PionZero())
      {
        if(++cpi0==randomPi0) 
	{
	  aCatch = true;
	  weight = npi0Weight;
	}
      }
      else
      {
        if(++cMeson==randomMeson) 
	{
	  aCatch = true;
	  weight = mesonWeight;
	}
      }
      if(aCatch)
      {
	buffer.push_back(*aSecTrack);
	buffer.back().SetWeight(aSecTrack->GetWeight()*weight);
      }
      else
      {
        delete aSecTrack;
      }
    }
    result->ClearSecondaries();
    // G4cerr << "pre"<<G4endl;
    result->AddSecondaries(buffer);
     // G4cerr << "bias exit"<<G4endl;
    
    return result;
  }
