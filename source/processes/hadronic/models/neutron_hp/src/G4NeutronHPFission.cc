// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPFission.hh"

  G4NeutronHPFission::G4NeutronHPFission()
  {
    SetMinEnergy( 0.0 );
    SetMaxEnergy( 20.*MeV );
    if(!getenv("NeutronHPCrossSections")) 
       G4Exception("Please setenv NeutronHPCrossSections to point to the neutron cross-section files.");
    dirName = getenv("NeutronHPCrossSections");
    G4String tString = "/Fission/";
    dirName = dirName + tString;
    numEle = G4Element::GetNumberOfElements();
    theFission = new G4NeutronHPChannel[numEle];

    for (G4int i=0; i<numEle; i++)
    { 
      if((*(G4Element::GetElementTable()))(i)->GetZ()>89)
      {
        theFission[i].Init((*(G4Element::GetElementTable()))(i), dirName);
        theFission[i].Register(&theFS);
      }
    }
  }
  
  G4NeutronHPFission::~G4NeutronHPFission()
  {
    delete [] theFission;
  }
  
  G4VParticleChange * G4NeutronHPFission::ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus)
  {
    G4Material * theMaterial = aTrack.GetMaterial();
    G4int n = theMaterial->GetNumberOfElements();
    xSec = new G4double[n];
    G4double sum=0;
    G4int i, it, index;
    const G4double * NumAtomsPerVolume = theMaterial->GetVecNbOfAtomsPerVolume();
    G4double rWeight;    
    for (i=0; i<n; i++)
    {
      index = theMaterial->GetElement(i)->GetIndex();
      rWeight = NumAtomsPerVolume[i];
      xSec[i] = theFission[index].GetXsec(aTrack.GetKineticEnergy());
      xSec[i] *= rWeight;
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
    return theFission[index].ApplyYourself(aTrack);
  }

