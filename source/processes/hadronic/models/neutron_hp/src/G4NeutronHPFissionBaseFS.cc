// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPFissionBaseFS.hh"
#include "G4ReactionProduct.hh"
#include "G4Nucleus.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4NeutronHPDataUsed.hh"

  G4NeutronHPFissionBaseFS::G4NeutronHPFissionBaseFS()
  { 
    hasXsec = true; 
    theXsection = new G4NeutronHPVector;
  }
  G4NeutronHPFissionBaseFS::~G4NeutronHPFissionBaseFS()
  {
    delete theXsection;
  }

  void G4NeutronHPFissionBaseFS::Init (G4double A, G4double Z, G4String & dirName, G4String & bit)
  {
    G4String tString = dirName;
    G4bool dbool;
    G4NeutronHPDataUsed aFile = theNames.GetName(A, Z, tString, bit, dbool);
    G4String filename = aFile.GetName();
    theBaseA = aFile.GetA();
    theBaseZ = aFile.GetZ();
    if(!dbool) 
    {
      hasAnyData = false;
      hasFSData = false; 
      hasXsec = false;
      return; // no data for exactly this isotope.
    }
    ifstream theData(filename, ios::in);
    G4int dummy;
    G4double dumm;
    if(!(theData))
    {
      hasFSData = false;
      hasXsec = false;
      hasAnyData = false;
      return; // no data for this FS for this isotope
    }
    theData >> dummy>>dummy;
    G4int total;
    theData >> total;
    theXsection->Init(theData, total, eV);
    if (!(theData >> dummy))
    {
      hasFSData = false;
      return;
    }
    theData >> dummy;

    theAngularDistribution.Init(theData);

    theData >> dummy >> dummy;

    theEnergyDistribution.Init(theData);
    
  }
  
G4DynamicParticleVector * G4NeutronHPFissionBaseFS::ApplyYourself(G4int nPrompt)
  {  
// if therere were no data for this isotope, break out.    
    if(!HasFSData()) return NULL;
    
    G4int i;
    G4DynamicParticleVector * aResult = new G4DynamicParticleVector;
    G4ReactionProduct boosted;
    boosted.Lorentz(theNeutron, theTarget);
    G4double eKinetic = boosted.GetKineticEnergy();
    
// Build neutrons
    G4ReactionProduct * theNeutrons = new G4ReactionProduct[nPrompt];
    for(i=0; i<nPrompt; i++)
    {
      theNeutrons[i].SetDefinition(G4Neutron::Neutron());
    }
    
// sample energies
   G4int dummy;
   for(i=0; i<nPrompt; i++)
   {
     // always in the lab system (if file-5)
     theNeutrons[i].SetKineticEnergy(theEnergyDistribution.Sample(eKinetic, dummy));
   }

// sample neutron angular distribution
   for(i=0; i<nPrompt; i++)
   {
     theAngularDistribution.SampleAndUpdate(theNeutrons[i]);
   }
   
// already in lab. Add neutrons to dynamic particle vector
   for(i=0; i<nPrompt; i++)
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
