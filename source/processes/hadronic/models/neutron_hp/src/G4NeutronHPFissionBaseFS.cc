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

  void G4NeutronHPFissionBaseFS::Init (G4double A, G4double Z, G4String & dirName, G4String & bit)
  {
    G4String tString = dirName;
    G4bool dbool;
    G4NeutronHPDataUsed aFile = theNames.GetName(A, Z, tString, bit, dbool);
    G4String filename = aFile.GetName();
    theBaseA = aFile.GetA();
    theBaseZ = aFile.GetZ();
    if(!dbool  || ( Z<2.5 && ( abs(theBaseZ - Z)>0.0001 || abs(theBaseA - A)>0.0001) ) )
    {
      hasAnyData = false;
      hasFSData = false; 
      hasXsec = false;
      return; // no data for exactly this isotope.
    }
#ifdef G4USE_STD_NAMESPACE
    G4std::ifstream theData(filename, G4std::ios::in);
#else
    ifstream theData(filename, ios::in|ios::nocreate);
#endif
    G4int dummy;
    G4double dumm;
    if(!(theData))
    {
      theData.close();
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
      theData.close();
      return;
    }
    theData >> dummy;

    theAngularDistribution.Init(theData);

    theData >> dummy >> dummy;

    theEnergyDistribution.Init(theData);
    theData.close();
    
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
