// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPInelasticCompFS.hh"
#include "G4Nucleus.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4Electron.hh"
#include "G4NeutronHPDataUsed.hh"
#include "G4ParticleTable.hh"

void G4NeutronHPInelasticCompFS::InitGammas(G4double AR, G4double ZR)
{
   char the[100] = {""};
   G4std::ostrstream ost(the, 100, G4std::ios::out);
   ost <<gammaPath<<"z"<<ZR<<".a"<<AR;
   G4String * aName = new G4String(the);
#ifdef G4USE_STD_NAMESPACE
   G4std::ifstream from(*aName, G4std::ios::in);
#else
   ifstream from(*aName, ios::in|ios::nocreate);
#endif
   if(!from) return; // no data found for this isotope
   G4std::ifstream theGammaData(*aName, G4std::ios::in);
    
   theGammas.Init(theGammaData);
   delete aName;
}

void G4NeutronHPInelasticCompFS::Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
{
  gammaPath = "/Inelastic/Gammas/";
    if(!getenv("NeutronHPCrossSections")) 
       G4Exception("Please setenv NeutronHPCrossSections to point to the neutron cross-section files.");
  G4String tBase = getenv("NeutronHPCrossSections");
  gammaPath = tBase+gammaPath;
  G4String tString = dirName;
  G4bool dbool;
  G4NeutronHPDataUsed aFile = theNames.GetName(A, Z, tString, aFSType, dbool);
  G4String filename = aFile.GetName();
    theBaseA = A;
    theBaseZ = G4int(Z+.5);
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
  if(!theData)
  {
    hasAnyData = false;
    hasFSData = false; 
    hasXsec = false;
    return;
  }
  // here we go
  G4int infoType, dataType, dummy;
  G4int sfType, it;
  hasFSData = false; 
  while (theData >> infoType)
  {
    hasFSData = true; 
    theData >> dataType;
    theData >> sfType >> dummy;
    it = 50;
    if(sfType>=600||(sfType<100&&sfType>=50)) it = sfType%50;
    if(dataType==3) 
    {
      theData >> dummy >> dummy;
      theXsection[it] = new G4NeutronHPVector;
      G4int total;
      theData >> total;
      theXsection[it]->Init(theData, total, eV);
    }
    else if(dataType==4)
    {
      theAngularDistribution[it] = new G4NeutronHPAngular;
      theAngularDistribution[it]->Init(theData);
    }
    else if(dataType==5)
    {
      theEnergyDistribution[it] = new G4NeutronHPEnergyDistribution;
      theEnergyDistribution[it]->Init(theData);
    }
    else if(dataType==6)
    {
      theEnergyAngData[it] = new G4NeutronHPEnAngCorrelation;
      theEnergyAngData[it]->Init(theData);
    }
    else if(dataType==12)
    {
      theFinalStatePhotons[it] = new G4NeutronHPPhotonDist;
      theFinalStatePhotons[it]->InitMean(theData);
    }
    else if(dataType==13)
    {
      theFinalStatePhotons[it] = new G4NeutronHPPhotonDist;
      theFinalStatePhotons[it]->InitPartials(theData);
    }
    else if(dataType==14)
    {
      theFinalStatePhotons[it]->InitAngular(theData);
    }
    else if(dataType==15)
    {
      theFinalStatePhotons[it]->InitEnergies(theData);
    }
    else
    {
      G4Exception("Data-type unknown to G4NeutronHPInelasticCompFS");
    }
  }
}

G4int G4NeutronHPInelasticCompFS::SelectExitChannel(G4double eKinetic)
{
  G4double running[50];
  running[0] = 0;
  G4int i;
  for(i=0; i<50; i++)
  {
    if(i!=0) running[i]=running[i-1];
    if(theXsection[i] != NULL) 
    {
      running[i] += G4std::max(0., theXsection[i]->GetXsec(eKinetic));
    }
  }
  G4double random = G4UniformRand();
  G4double sum = running[49];
  G4int it = 50;
  if(0!=sum)
  {
    for(i=0; i<50; i++)
    {
      it = i;
      if(random < running[i]/sum) break;
    }
  }
//debug:  it = 1;
  return it;
}

void G4NeutronHPInelasticCompFS::CompositeApply(const G4Track & theTrack, G4ParticleDefinition * aDefinition)
{
    theResult.Initialize(theTrack); 

// prepare neutron
    G4double eKinetic = theTrack.GetKineticEnergy();
    const G4DynamicParticle *incidentParticle = theTrack.GetDynamicParticle();
    G4ReactionProduct theNeutron( incidentParticle->GetDefinition() );
    theNeutron.SetMomentum( incidentParticle->GetMomentum() );
    theNeutron.SetKineticEnergy( eKinetic );

// prepare target
    G4int i;
    for(i=0; i<50; i++) if(theXsection[i] != NULL) break; 
    G4double targetMass=0;
    G4double eps = 0.0001;
    targetMass = ( G4NucleiPropertiesTable::GetNuclearMass(theBaseZ+eps, theBaseA+eps)) /
                   G4Neutron::Neutron()->GetPDGMass();
//    if(theEnergyAngData[i]!=NULL)
//        targetMass = theEnergyAngData[i]->GetTargetMass();
//    else if(theAngularDistribution[i]!=NULL)
//        targetMass = theAngularDistribution[i]->GetTargetMass();
//    else if(theFinalStatePhotons[50]!=NULL)
//        targetMass = theFinalStatePhotons[50]->GetTargetMass();
    G4Nucleus aNucleus;
    G4ReactionProduct theTarget; 
    theTarget = aNucleus.GetThermalNucleus(targetMass);

// prepare the residual mass
    G4double residualMass=0;
    G4double residualZ = theBaseZ - aDefinition->GetPDGCharge();
    G4double residualA = theBaseA - aDefinition->GetBaryonNumber()+1;
    residualMass = ( G4NucleiPropertiesTable::GetNuclearMass(residualZ+eps, residualA+eps) ) /
                     G4Neutron::Neutron()->GetPDGMass();

// prepare energy in target rest frame
    G4ReactionProduct boosted;
    boosted.Lorentz(theNeutron, theTarget);
    eKinetic = boosted.GetKineticEnergy();
    G4double momentumInCMS = boosted.GetTotalMomentum();
  
// select exit channel for composite FS class.
    G4int it = SelectExitChannel(eKinetic);
   
// set target and neutron in the relevant exit channel
    InitDistributionInitialState(theNeutron, theTarget, it);    

    G4ReactionProductVector * thePhotons = NULL;
    G4ReactionProductVector * theParticles = NULL;
    G4ReactionProduct aHadron;
    aHadron.SetDefinition(aDefinition); // what if only cross-sections exist ==> Na 23 11 @@@@    
    G4double availableEnergy = theNeutron.GetKineticEnergy() + theNeutron.GetMass() - aHadron.GetMass() +
                             (targetMass - residualMass)*G4Neutron::Neutron()->GetPDGMass();
    G4int nothingWasKnownOnHadron = 0;
    G4int dummy;
    G4int nSecGamma = 0;
    G4double eGamm = 0;
    G4int iLevel=it-1;
    if(50==it) 
    {
      iLevel=-1;
      aHadron.SetKineticEnergy(availableEnergy*residualMass*G4Neutron::Neutron()->GetPDGMass()/
                               (aHadron.GetMass()+residualMass*G4Neutron::Neutron()->GetPDGMass()));
      aHadron.SetMomentum(theNeutron.GetMomentum()*(1./theNeutron.GetTotalMomentum())*
                           sqrt(aHadron.GetTotalEnergy()*aHadron.GetTotalEnergy()-
                                aHadron.GetMass()*aHadron.GetMass()));
    }
    else
    {
      while( iLevel!=-1 && theGammas.GetLevel(iLevel)==NULL ) iLevel--;
    }
    if(theAngularDistribution[it]!= NULL)
    {
      if(theEnergyDistribution[it]!=NULL)
      {
        aHadron.SetKineticEnergy(theEnergyDistribution[it]->Sample(eKinetic, dummy));
        G4double eSecN = aHadron.GetKineticEnergy();
        eGamm = eKinetic-eSecN;
        for(iLevel=theGammas.GetNumberOfLevels()-1; iLevel>=0; iLevel--)
        {
          if(theGammas.GetLevelEnergy(iLevel)<eGamm) break;
        }
        G4double random = 2*G4UniformRand();
        iLevel+=G4int(random);
        if(iLevel>theGammas.GetNumberOfLevels()-1)iLevel = theGammas.GetNumberOfLevels()-1;
      }
      else
      {
        G4double eExcitation = 0;
        if(iLevel>=0) eExcitation = theGammas.GetLevel(iLevel)->GetLevelEnergy();    
        
	aHadron.SetKineticEnergy(aHadron.GetKineticEnergy() - eExcitation);
	// consistency of data assumed....@@@@@
	
      }
      theAngularDistribution[it]->SampleAndUpdate(aHadron);
      if(theFinalStatePhotons[it] == NULL)
      {
	thePhotons = theGammas.GetDecayGammas(iLevel);
	eGamm -= theGammas.GetLevelEnergy(iLevel);
	if(eGamm>0) // @ ok for now, but really needs an efficient way of correllated sampling @
	{
          G4ReactionProduct * theRestEnergy = new G4ReactionProduct;
          theRestEnergy->SetDefinition(G4Gamma::Gamma());
          theRestEnergy->SetKineticEnergy(eGamm);
          G4double costh = 2.*G4UniformRand()-1.;
          G4double phi = twopi*G4UniformRand();
          theRestEnergy->SetMomentum(eGamm*sin(acos(costh))*cos(phi), 
                                     eGamm*sin(acos(costh))*sin(phi),
                                     eGamm*costh);
          if(thePhotons == NULL) thePhotons = new G4ReactionProductVector;
          thePhotons->insert(theRestEnergy);
	}
      }
    }
    else if(theEnergyAngData[it] != NULL)  
    {
      theParticles = theEnergyAngData[it]->Sample(eKinetic);
    }
    else
    {
      // @@@ what to do, if we have photon data, but no info on the hadron itself
      nothingWasKnownOnHadron = 1;
    }
    if(theFinalStatePhotons[it]!=NULL) 
    {
      // the photon distributions are in the Nucleus rest frame.
      G4ReactionProduct boosted;
      boosted.Lorentz(theNeutron, theTarget);
      G4double anEnergy = boosted.GetKineticEnergy();
      thePhotons = theFinalStatePhotons[it]->GetPhotons(anEnergy);
      G4double aBaseEnergy = theFinalStatePhotons[it]->GetLevelEnergy();
      G4double testEnergy = 0;
      if(thePhotons!=NULL && thePhotons->entries()!=0) aBaseEnergy-=thePhotons->at(0)->GetTotalEnergy();
      if(theFinalStatePhotons[it]->NeedsCascade())
      {
	while(aBaseEnergy>0.01*keV)
        {
          // cascade down the levels
	  G4bool foundMatchingLevel = false;
          G4int closest;
	  G4double deltaEold = -1;
	  for(G4int i=1; i<it; i++)
          {
            if(theFinalStatePhotons[i]!=NULL) 
            {
              testEnergy = theFinalStatePhotons[i]->GetLevelEnergy();
            }
            else
            {
              testEnergy = 0;
            }
	    G4double deltaE = abs(testEnergy-aBaseEnergy);
            if(deltaE<0.1*keV)
            {
              G4ReactionProductVector * theNext = 
        	theFinalStatePhotons[i]->GetPhotons(anEnergy);
              thePhotons->insert(theNext->at(0));
              aBaseEnergy = testEnergy-theNext->at(0)->GetTotalEnergy();
              delete theNext;
	      foundMatchingLevel = true;
              break; // ===>
            }
	    if(deltaE<deltaEold||deltaEold<0.)
	    {
	      closest = i;
	      deltaEold = deltaE;     
	    }
          } // <=== the break goes here.
	  if(!foundMatchingLevel)
	  {
            G4ReactionProductVector * theNext = 
               theFinalStatePhotons[closest]->GetPhotons(anEnergy);
            thePhotons->insert(theNext->at(0));
	    testEnergy = theFinalStatePhotons[closest]->GetLevelEnergy();
            aBaseEnergy = testEnergy-theNext->at(0)->GetTotalEnergy();
            delete theNext;
	  }
        } 
      }
    }
    if(thePhotons!=NULL)
    {
      for(i=0; i<thePhotons->length(); i++)
      {
	// back to lab
	thePhotons->at(i)->Lorentz(*(thePhotons->at(i)), -1.*theTarget);
      }
    }
    if(nothingWasKnownOnHadron)
    {
      G4double totalPhotonEnergy = 0;
      if(thePhotons!=NULL)
      {
        G4int nPhotons = thePhotons->length();
	for(i=0; i<nPhotons; i++)
        {
          totalPhotonEnergy += thePhotons->at(i)->GetTotalEnergy();
        }
      }
      availableEnergy -= totalPhotonEnergy;
      residualMass += totalPhotonEnergy/G4Neutron::Neutron()->GetPDGMass();
      aHadron.SetKineticEnergy(availableEnergy*residualMass*G4Neutron::Neutron()->GetPDGMass()/
                               (aHadron.GetMass()+residualMass*G4Neutron::Neutron()->GetPDGMass()));
      G4double CosTheta = 1.0 - 2.0*G4UniformRand();
      G4double SinTheta = sqrt(1.0 - CosTheta*CosTheta);
      G4double Phi = twopi*G4UniformRand();
      G4ThreeVector Vector(cos(Phi)*SinTheta, sin(Phi)*SinTheta, CosTheta);
      aHadron.SetMomentum(Vector* sqrt(aHadron.GetTotalEnergy()*aHadron.GetTotalEnergy()-
                                       aHadron.GetMass()*aHadron.GetMass()));
    }

// fill the result
// Beware - the recoil is not necessarily in the particles...
// Can be calculated from momentum conservation?
// The idea is that the particles ar emitted forst, and the gammas only once the
// recoil is on the residual; assumption is that gammas do not contribute to 
// the recoil.
// This needs more design @@@

    G4int nSecondaries = 2; // the hadron and the recoil
    G4bool needsSeparateRecoil = false;
    G4int totalBaryonNumber = 0;
    G4int totalCharge = 0;
    G4ThreeVector totalMomentum(0);
    if(theParticles != NULL) 
    {
      nSecondaries = theParticles->length();
      G4ParticleDefinition * aDef;
      for(i=0; i<theParticles->length(); i++)
      {
        aDef = theParticles->at(i)->GetDefinition();
	totalBaryonNumber+=aDef->GetBaryonNumber();
	totalCharge+=G4int(aDef->GetPDGCharge()+eps);
        totalMomentum += theParticles->at(i)->GetMomentum();
      } 
      if(totalBaryonNumber!=G4int(theBaseA+eps+incidentParticle->GetDefinition()->GetBaryonNumber())) 
      {
        needsSeparateRecoil = true;
	nSecondaries++;
	residualA = G4int(theBaseA+eps+incidentParticle->GetDefinition()->GetBaryonNumber()
	                  -totalBaryonNumber);
	residualZ = G4int(theBaseZ+eps+incidentParticle->GetDefinition()->GetPDGCharge()
	                  -totalCharge);
      }
    }
    
    G4int nPhotons = 0;
    if(thePhotons!=NULL) nPhotons = thePhotons->length();
    nSecondaries += nPhotons;
    
    theResult.SetNumberOfSecondaries(nSecondaries);
    
    G4DynamicParticle * theSec;
    
    if( theParticles==NULL )
    {
      theSec = new G4DynamicParticle;   
      theSec->SetDefinition(aHadron.GetDefinition());
      theSec->SetMomentum(aHadron.GetMomentum());
      theResult.AddSecondary(theSec);    
 
 	aHadron.Lorentz(aHadron, theTarget);
        G4ReactionProduct theResidual;   
        theResidual.SetDefinition(G4ParticleTable::GetParticleTable()->GetIon(residualZ, residualA, 0));  
        theResidual.SetKineticEnergy(aHadron.GetKineticEnergy()*aHadron.GetMass()/theResidual.GetMass());
        theResidual.SetMomentum(-1.*aHadron.GetMomentum());
	theResidual.Lorentz(theResidual, -1.*theTarget);
	G4ThreeVector totalPhotonMomentum(0,0,0);
	if(thePhotons!=NULL)
	{
          for(i=0; i<nPhotons; i++)
          {
            totalPhotonMomentum += thePhotons->at(i)->GetMomentum();
          }
	}
        theSec = new G4DynamicParticle;   
        theSec->SetDefinition(theResidual.GetDefinition());
        theSec->SetMomentum(theResidual.GetMomentum()-totalPhotonMomentum);
        theResult.AddSecondary(theSec);    
    }
    else
    {
      for(i=0; i<theParticles->length(); i++)
      {
        theSec = new G4DynamicParticle; 
        theSec->SetDefinition(theParticles->at(i)->GetDefinition());
        theSec->SetMomentum(theParticles->at(i)->GetMomentum());
        theResult.AddSecondary(theSec); 
        delete theParticles->at(i); 
      } 
      delete theParticles;
      if(needsSeparateRecoil)
      {
        G4ReactionProduct theResidual;   
        theResidual.SetDefinition(G4ParticleTable::GetParticleTable()->GetIon(residualZ, residualA, 0));  
        G4double resiualKineticEnergy  = theResidual.GetMass()*theResidual.GetMass();
                 resiualKineticEnergy += totalMomentum*totalMomentum;
  	         resiualKineticEnergy  = sqrt(resiualKineticEnergy) - theResidual.GetMass();
//        cout << "Kinetic energy of the residual = "<<resiualKineticEnergy<<endl;
	theResidual.SetKineticEnergy(resiualKineticEnergy);
        theResidual.SetMomentum(-1.*totalMomentum);
        theSec = new G4DynamicParticle;   
        theSec->SetDefinition(theResidual.GetDefinition());
        theSec->SetMomentum(theResidual.GetMomentum());
        theResult.AddSecondary(theSec);  
      }  
    }
    if(thePhotons!=NULL)
    {
      for(i=0; i<nPhotons; i++)
      {
        theSec = new G4DynamicParticle;    
        theSec->SetDefinition(G4Gamma::Gamma());
        theSec->SetMomentum(thePhotons->at(i)->GetMomentum());
        theResult.AddSecondary(theSec); 
        delete thePhotons->at(i);
      }
// some garbage collection
      delete thePhotons;
    }
// clean up the primary neutron
    theResult.SetStatusChange(fStopAndKill);
}
