// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureFS.hh"
#include "G4NeutronHPDeExGammas.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

  G4NeutronHPCapture::G4NeutronHPCapture()
  {
//    G4cout << "Capture : start of construction!!!!!!!!"<<G4endl;
    if(!getenv("NeutronHPCrossSections")) 
       G4Exception("Please setenv NeutronHPCrossSections to point to the neutron cross-section files.");
    dirName = getenv("NeutronHPCrossSections");
    G4String tString = "/Capture/";
    dirName = dirName + tString;
    numEle = G4Element::GetNumberOfElements();
//    G4cout << "+++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl;
//    G4cout <<"Disname="<<dirName<<" numEle="<<numEle<<G4endl;
    theCapture = new G4NeutronHPChannel[numEle];
//    G4cout <<"G4NeutronHPChannel constructed"<<G4endl;
    G4NeutronHPCaptureFS * theFS = new G4NeutronHPCaptureFS;
    for (G4int i=0; i<numEle; i++)
    {
//      G4cout << "initializing theCapture "<<i<<" "<< numEle<<G4endl;
      theCapture[i].Init((*(G4Element::GetElementTable()))(i), dirName);
      theCapture[i].Register(theFS);
    }
    delete theFS;
//    G4cout << "-------------------------------------------------"<<G4endl;
//    G4cout << "Leaving G4NeutronHPCapture::G4NeutronHPCapture"<<G4endl;
  }
  
  G4NeutronHPCapture::~G4NeutronHPCapture()
  {
    delete [] theCapture;
//    G4cout << "Leaving G4NeutronHPCapture::~G4NeutronHPCapture"<<G4endl;
  }
  
  G4VParticleChange * G4NeutronHPCapture::ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus)
  {
    G4Material * theMaterial = aTrack.GetMaterial();
    G4int n = theMaterial->GetNumberOfElements();
    xSec = new G4double[n];
    G4double sum=0;
    G4int i, index;
    const G4double * NumAtomsPerVolume = theMaterial->GetVecNbOfAtomsPerVolume();
    G4double rWeight;    
    for (i=0; i<n; i++)
    {
      index = theMaterial->GetElement(i)->GetIndex();
      rWeight = NumAtomsPerVolume[i];
      xSec[i] = theCapture[index].GetXsec(aTrack.GetKineticEnergy());
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
    if(i==n) i=G4std::max(0, n-1);
    delete [] xSec;
    return theCapture[index].ApplyYourself(aTrack);
  }
