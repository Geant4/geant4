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

#include "G4HadLeadBias.hh"
#include "G4Gamma.hh"
#include "G4PionZero.hh"
#include "Randomize.hh"
#include "G4ParticleChange.hh"

  G4VParticleChange * G4HadLeadBias::Bias(G4VParticleChange * result)
  {
    G4cerr << "bias enter"<<G4endl;
    G4int nMeson(0), nBaryon(0), npi0(0), ngamma(0), nLepton(0);
    G4int i(0);
    G4int maxE = -1;
    G4double emax = 0;
    G4ParticleChange * temp;
    if(result->GetStatusChange()==fAlive) 
    {
      temp = dynamic_cast<G4ParticleChange *>(result);
      if(temp) emax = temp->GetEnergyChange();
    }
    G4cout << "max energy "<<G4endl;
    for(i=0;i<result->GetNumberOfSecondaries();i++)
    {
      if(result->GetSecondary(i)->GetKineticEnergy()>emax)
      {
        maxE = i;
	emax = result->GetSecondary(i)->GetKineticEnergy();
      }
    }
    G4cout <<"loop1"<<G4endl;
    for(i=0; i<result->GetNumberOfSecondaries(); i++)
    {
      G4Track* aSecTrack = result->GetSecondary(i);
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
     G4cout << "BiasDebug 1 = "<<result->GetNumberOfSecondaries()<<" "
            <<nMeson<<" "<< nBaryon<<" "<< npi0<<" "<< ngamma<<" "<< nLepton<<G4endl;
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
    
    G4std::vector<G4Track*> buffer;
    G4int cMeson(0), cBaryon(0), cpi0(0), cgamma(0), cLepton(0);
    for(i=0; i<result->GetNumberOfSecondaries(); i++)
    {
      G4bool aCatch = false;
      G4double weight = 1;
      G4Track* aSecTrack = result->GetSecondary(i);
      if(i==maxE)
      {
        aCatch = true;
	weight = 1;
      }
      else if(aSecTrack->GetDefinition()->GetBaryonNumber()!=0) 
      {
	if(++cBaryon==randomBaryon) 
	{
	  aCatch = true;
	  weight = baryonWeight;
	}
      }
      else if(aSecTrack->GetDefinition()->GetLeptonNumber()!=0) 
      {
        if(++cLepton==randomLepton) 
	{
	  aCatch = true;
	  weight = leptonWeight;
	}
      }
      else if(aSecTrack->GetDefinition()==G4Gamma::Gamma())
      {
        if(++cgamma==randomGamma) 
	{
	  aCatch = true;
	  weight = gammaWeight;
	}
      }
      else if(aSecTrack->GetDefinition()==G4PionZero::PionZero())
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
	buffer.push_back(aSecTrack);
	aSecTrack->SetWeight(aSecTrack->GetWeight()*weight);
      }
      else
      {
        delete aSecTrack;
      }
    }
    result->Clear();
    result->SetNumberOfSecondaries(buffer.size());
    // G4cerr << "pre"<<G4endl;
    for(i=0;i<static_cast<G4int>(buffer.size());i++)
    {
      result->AddSecondary(buffer[i]);
    }
     G4cerr << "bias exit"<<G4endl;
    
    return result;
  }
