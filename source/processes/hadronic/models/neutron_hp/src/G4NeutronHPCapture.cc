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
    delete [] xSec;
    if(aTrack.GetKineticEnergy()<100*keV)
    {
      G4NeutronHPDeExGammas theGammas;
      G4int aA = theMaterial->GetElement(i)->GetN();
      G4int aZ = theMaterial->GetElement(i)->GetZ();
      char the[100] = {""};
      G4std::ostrstream ost(the, 100, G4std::ios::out);
      ost
      <<getenv("NeutronHPCrossSections")<<"/Inelastic/Gammas/"<<"z"<<aZ<<".a"<<aA;
      G4String * aName = new G4String(the);
#ifdef G4USE_STD_NAMESPACE
      G4std::ifstream from(*aName, G4std::ios::in);
#else
      ifstream from(*aName, ios::in|ios::nocreate);
#endif
      G4std::ifstream theGammaData(*aName, G4std::ios::in);
     
      theGammas.Init(theGammaData);
      G4double theGammaEnergy = aTrack.GetKineticEnergy();
      theGammaEnergy+=G4Neutron::Neutron()->GetPDGMass();
      theGammaEnergy+=
             G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(aZ,aA);
      theGammaEnergy-=
             G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(aZ,aA+1);
    
      if(theGammaEnergy<1.1*theGammas.GetLevelEnergy(theGammas.GetNumberOfLevels()-1))
      {
        G4ReactionProductVector * thePhotons = 0;
        G4ReactionProductVector * theOtherPhotons = NULL;
        G4int iLevel;
        while(theGammaEnergy>=theGammas.GetLevelEnergy(0))
        {
          for(iLevel=theGammas.GetNumberOfLevels()-1; iLevel>=0; iLevel--)
          {
    	    if(theGammas.GetLevelEnergy(iLevel)<theGammaEnergy) break;
          }
          if(iLevel==0||iLevel==theGammas.GetNumberOfLevels()-1)
          {
 	    theOtherPhotons = theGammas.GetDecayGammas(iLevel);
          }
          else
          {
  	    G4double random = G4UniformRand();
	    G4double eLow  = theGammas.GetLevelEnergy(iLevel);
	    G4double eHigh = theGammas.GetLevelEnergy(iLevel+1);
	    if(random > (eHigh-eLow)/(theGammaEnergy-eLow)) iLevel++;
	    theOtherPhotons = theGammas.GetDecayGammas(iLevel);
          }
          if(thePhotons==NULL) thePhotons = new G4ReactionProductVector;
          if(theOtherPhotons != NULL)
          {
            for(G4int ii=0; ii<theOtherPhotons->length(); ii++)
            {
              thePhotons->insert(theOtherPhotons->at(ii));
            }
            delete theOtherPhotons; 
          }
          theGammaEnergy -= theGammas.GetLevelEnergy(iLevel);
          if(iLevel == -1) break;
        }  
      
        // clean up the primary neutron
        theResult.SetStatusChange(fStopAndKill);
      
        // fill particle change
        theResult.Initialize(aTrack);
        G4int nSecondaries = thePhotons->length();
        theResult.SetNumberOfSecondaries(nSecondaries);
        G4DynamicParticle * theSec;
        for(G4int gammaCount=0; gammaCount<nSecondaries; gammaCount++)
        {
          theSec = new G4DynamicParticle;    
          theSec->SetDefinition(thePhotons->at(gammaCount)->GetDefinition());
          theSec->SetMomentum(thePhotons->at(gammaCount)->GetMomentum());
          theResult.AddSecondary(theSec); 
          delete thePhotons->at(gammaCount);
        }
        delete thePhotons;
      
        // return   
        if(0!=nSecondaries) return &theResult;
      }
      return theCapture[index].ApplyYourself(aTrack);
    }
    else
    {
      return theCapture[index].ApplyYourself(aTrack);
    }
  }
