// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticFS.hh"

  G4NeutronHPElastic::G4NeutronHPElastic()
  {
    G4NeutronHPElasticFS * theFS = new G4NeutronHPElasticFS;
    if(!getenv("NeutronHPCrossSections")) 
       G4Exception("Please setenv NeutronHPCrossSections to point to the neutron cross-section files.");
    dirName = getenv("NeutronHPCrossSections");
    G4String tString = "/Elastic/";
    dirName = dirName + tString;
//    G4cout <<"G4NeutronHPElastic::G4NeutronHPElastic testit "<<dirName<<endl;
    numEle = G4Element::GetNumberOfElements();
    theElastic = new G4NeutronHPChannel[numEle];
    for (G4int i=0; i<numEle; i++)
    {
      theElastic[i].Init((*(G4Element::GetElementTable()))(i), dirName);
      while(!theElastic[i].Register(theFS));
    }
    delete theFS;
  }
  
  G4NeutronHPElastic::~G4NeutronHPElastic()
  {
    delete [] theElastic;
  }
  
  G4VParticleChange * G4NeutronHPElastic::ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus)
  {
    G4Material * theMaterial = aTrack.GetMaterial();
    G4int n = theMaterial->GetNumberOfElements();
    G4int it = n-1;
    xSec = new G4double[n];
    G4double sum=0;
    G4int i, index;
    for (i=0; i<n; i++)
    {
      index = theMaterial->GetElement(i)->GetIndex();
      xSec[i] = theElastic[index].GetXsec(aTrack.GetKineticEnergy());
      sum+=xSec[i];
    }
    G4double random = G4UniformRand();
    G4double running = 0;
    for (i=0; i<n; i++)
    {
      running += xSec[i];
      index = theMaterial->GetElement(i)->GetIndex();
      if(random<=running/sum) break;
    }
    delete [] xSec;
     // it is element-wise initialised.
    return theElastic[index].ApplyYourself(aTrack); 
  }
