// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPFSFissionFS.hh"
#include "G4ReactionProduct.hh"
#include "G4Nucleus.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4NeutronHPDataUsed.hh"

  void G4NeutronHPFSFissionFS::Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
  {
    G4String tString = "/FS/";
    G4bool dbool;
    G4NeutronHPDataUsed aFile = theNames.GetName(A, Z, dirName, tString, dbool);
    G4String filename = aFile.GetName();
    if(!dbool)
    {
      hasAnyData = false;
      hasFSData = false; 
      hasXsec = false;
      return;
    }
#ifdef G4USE_STD_NAMESPACE
    G4std::ifstream theData(filename, G4std::ios::in);
#else
    ifstream theData(filename, ios::in|ios::nocreate);
#endif
    // here it comes
    G4int infoType, dataType;
    hasFSData = false; 
    while (theData >> infoType)
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
          G4cout << "G4NeutronHPFSFissionFS::Init: unknown data type"<<dataType<<G4endl;
          G4Exception("G4NeutronHPFSFissionFS::Init: unknown data type");
          break;
      }
    }
    targetMass = theFinalStateNeutrons.GetTargetMass();
  }
  
  
  G4DynamicParticleVector * G4NeutronHPFSFissionFS::ApplyYourself(G4int nPrompt, 
                                                 G4int nDelayed, G4double * theDecayConst)
  {  
    G4int i;
    G4DynamicParticleVector * aResult = new G4DynamicParticleVector;
    G4ReactionProduct boosted;
    boosted.Lorentz(theNeutron, theTarget);
    G4double eKinetic = boosted.GetKineticEnergy();
    
// Build neutrons
    G4ReactionProduct * theNeutrons = new G4ReactionProduct[nPrompt+nDelayed];
    for(i=0; i<nPrompt+nDelayed; i++)
    {
      theNeutrons[i].SetDefinition(G4Neutron::Neutron());
    }
    
// sample energies
   G4int it, dummy;
   G4double tempE;
   for(i=0; i<nPrompt; i++)
   {
     tempE = thePromptNeutronEnDis.Sample(eKinetic, dummy); // energy distribution (file5) always in lab
     theNeutrons[i].SetKineticEnergy(tempE);
   }
   for(i=nPrompt; i<nPrompt+nDelayed; i++)
   {
     theNeutrons[i].SetKineticEnergy(theDelayedNeutronEnDis.Sample(eKinetic, it));  // dito
     if(it==0) theNeutrons[i].SetKineticEnergy(thePromptNeutronEnDis.Sample(eKinetic, dummy));
     theDecayConst[i-nPrompt] = theFinalStateNeutrons.GetDecayConstant(it); // this is returned
   }

// sample neutron angular distribution
   for(i=0; i<nPrompt+nDelayed; i++)
   {
     theNeutronAngularDis.SampleAndUpdate(theNeutrons[i]); // angular comes back in lab automatically
   }
   
// already in lab. Add neutrons to dynamic particle vector
   for(i=0; i<nPrompt+nDelayed; i++)
   {
      G4DynamicParticle * it = new G4DynamicParticle;
      it->SetDefinition(theNeutrons[i].GetDefinition());
      it->SetMomentum(theNeutrons[i].GetMomentum());
      aResult->insert(it);
   }
   delete [] theNeutrons;
// return the result
   return aResult;
  }

void G4NeutronHPFSFissionFS::SampleNeutronMult(G4int&all, G4int&Prompt, G4int&delayed, G4double eKinetic, G4int off)
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
     all = RandPoisson::shoot(totalNeutronMulti-off);
     all += off;
   }
   else
   {   
     Prompt  = RandPoisson::shoot(promptNeutronMulti-off);
     Prompt += off;
     delayed = RandPoisson::shoot(delayedNeutronMulti);
     all = Prompt+delayed;
   }
}

G4DynamicParticleVector * G4NeutronHPFSFissionFS::GetPhotons()
{
// sample photons
   G4ReactionProductVector * temp;
   G4ReactionProduct boosted;
// the photon distributions are in the Nucleus rest frame.
   boosted.Lorentz(theNeutron, theTarget);
   G4double anEnergy = boosted.GetKineticEnergy();
   temp = theFinalStatePhotons.GetPhotons(anEnergy);
   if(temp == NULL) return NULL;

// lorentz transform, and add photons to final state
   G4int i;
   G4DynamicParticleVector * result = new G4DynamicParticleVector;
   for(i=0; i<temp->length(); i++)
   {
     // back to lab
     temp->at(i)->Lorentz(*temp->at(i), -1.*theTarget);
     G4DynamicParticle * theOne = new G4DynamicParticle;
     theOne->SetDefinition(temp->at(i)->GetDefinition());
     theOne->SetMomentum(temp->at(i)->GetMomentum());
     result->insert(theOne);
     delete temp->at(i);
   }
   delete temp;
   return result;
}
