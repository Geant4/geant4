// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPInelasticBaseFS.hh"
#include "G4Nucleus.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4Electron.hh"
#include "G4NeutronHPDataUsed.hh"

   G4NeutronHPInelasticBaseFS::G4NeutronHPInelasticBaseFS()
  {
    hasXsec = true; 
    theXsection = new G4NeutronHPVector;
    
    theEnergyDistribution = NULL;
    theFinalStatePhotons = NULL;
    theEnergyAngData = NULL;
    theAngularDistribution = NULL;
  }
   G4NeutronHPInelasticBaseFS::~G4NeutronHPInelasticBaseFS()
  {
    delete theXsection;
    if(theEnergyDistribution!=NULL) delete theEnergyDistribution;
    if(theFinalStatePhotons!=NULL) delete theFinalStatePhotons;
    if(theEnergyAngData!=NULL) delete theEnergyAngData;
    if(theAngularDistribution!=NULL) delete theAngularDistribution;
  }

void G4NeutronHPInelasticBaseFS::InitGammas(G4double AR, G4double ZR)
{
   char the[100] = {""};
   ostrstream ost(the, 100, ios::out);
   ost <<gammaPath<<"z"<<ZR<<".a"<<AR;
   G4String * aName = new G4String(the);
   ifstream from(*aName, ios::in);
   if(!from) return; // no data found for this isotope
   ifstream theGammaData(*aName, ios::in);
    
   G4double eps = 0.001;
   theNuclearMassDifference = 
       G4NucleiPropertiesTable::GetBindingEnergy(ZR+eps,AR+eps) -
       G4NucleiPropertiesTable::GetBindingEnergy(theBaseZ+eps, theBaseA+eps);
   theGammas.Init(theGammaData);
   delete aName;
}

void G4NeutronHPInelasticBaseFS::Init (G4double A, G4double Z, G4String & dirName, G4String & bit)
{
  gammaPath = "/Inelastic/Gammas/";
  G4String tBase = getenv("NeutronHPCrossSections");
  gammaPath = tBase+gammaPath;
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
    return;
  }
  ifstream theData(filename, ios::in);
  if(!(theData))
  {
    hasAnyData = false;
    hasFSData = false; 
    hasXsec = false;
    return; // no data for exactly this isotope and FS
  }
  // here we go
  G4int infoType, dataType, dummy=INT_MAX;
  hasFSData = false; 
  while (theData >> infoType)
  {
    theData >> dataType;
    if(dummy==INT_MAX) theData >> dummy >> dummy;
    if(dataType==3) 
    {
      G4int total;
      theData >> total;
      theXsection->Init(theData, total, eV);
    }
    else if(dataType==4)
    {
      theAngularDistribution = new G4NeutronHPAngular;
      theAngularDistribution->Init(theData);
      hasFSData = true; 
    }
    else if(dataType==5)
    {
      theEnergyDistribution = new G4NeutronHPEnergyDistribution;
      theEnergyDistribution->Init(theData);
      hasFSData = true; 
    }
    else if(dataType==6)
    {
      theEnergyAngData = new G4NeutronHPEnAngCorrelation;
      theEnergyAngData->Init(theData);
      hasFSData = true; 
    }
    else if(dataType==12)
    {
      theFinalStatePhotons = new G4NeutronHPPhotonDist;
      theFinalStatePhotons->InitMean(theData);
      hasFSData = true; 
    }
    else if(dataType==13)
    {
      theFinalStatePhotons = new G4NeutronHPPhotonDist;
      theFinalStatePhotons->InitPartials(theData);
      hasFSData = true; 
    }
    else if(dataType==14)
    {
      theFinalStatePhotons->InitAngular(theData);
      hasFSData = true; 
    }
    else if(dataType==15)
    {
      theFinalStatePhotons->InitEnergies(theData);
      hasFSData = true; 
    }
    else
    {
      G4Exception("Data-type unknown to G4NeutronHPInelasticBaseFS");
    }
  }
}
  
void G4NeutronHPInelasticBaseFS::BaseApply(const G4Track & theTrack, 
                                           G4ParticleDefinition ** theDefs, 
                                           G4int nDef)
{
  theResult.Initialize(theTrack); 

// prepare neutron
  G4double eKinetic = theTrack.GetKineticEnergy();
  const G4DynamicParticle *incidentParticle = theTrack.GetDynamicParticle();
  G4ReactionProduct theNeutron( incidentParticle->GetDefinition() );
  theNeutron.SetMomentum( incidentParticle->GetMomentum() );
  theNeutron.SetKineticEnergy( eKinetic );

// prepare target
  G4double targetMass;
  G4double eps = 0.0001;
  targetMass = ( G4NucleiPropertiesTable::GetAtomicMass(theBaseZ+eps, theBaseA+eps)-
                            theBaseZ*G4Electron::ElectronDefinition()->GetPDGMass() ) /
               G4Neutron::Neutron()->GetPDGMass();
  if(theEnergyAngData!=NULL)
      targetMass = theEnergyAngData->GetTargetMass();
  if(theAngularDistribution!=NULL)
      targetMass = theAngularDistribution->GetTargetMass();
  G4Nucleus aNucleus;
  G4ReactionProduct theTarget; 
  theTarget = aNucleus.GetThermalNucleus(targetMass);

// prepare energy in target rest frame
  G4ReactionProduct boosted;
  boosted.Lorentz(theNeutron, theTarget);
  eKinetic = boosted.GetKineticEnergy();
  
// Take N-body phase-space distribution, if no other data present.
  if(!HasFSData()) // adding the residual is trivial here @@@
  {
    G4NeutronHPNBodyPhaseSpace thePhaseSpaceDistribution;
    G4double aPhaseMass=0;
    G4int ii;
    for(ii=0; ii<nDef; ii++) 
    {
      aPhaseMass+=theDefs[ii]->GetPDGMass();
    }
    theResult.SetNumberOfSecondaries(nDef);
    thePhaseSpaceDistribution.Init(aPhaseMass, nDef);
    thePhaseSpaceDistribution.SetNeutron(&theNeutron);
    thePhaseSpaceDistribution.SetTarget(&theTarget);
    for(ii=0; ii<nDef; ii++) 
    {
      G4double massCode = 1000.*abs(theDefs[ii]->GetPDGCharge());
      massCode += theDefs[ii]->GetBaryonNumber(); 
      G4double dummy = 0;
      G4ReactionProduct * aSec = thePhaseSpaceDistribution.Sample(eKinetic, massCode, dummy);
      aSec->Lorentz(*aSec, -1.*theTarget);
      G4DynamicParticle * aPart = new G4DynamicParticle();
      aPart->SetDefinition(aSec->GetDefinition());
      aPart->SetMomentum(aSec->GetMomentum());
      delete aSec;
      theResult.AddSecondary(aPart);     
    }   
    theResult.SetStatusChange(fStopAndKill);
    return;
  }

// set target and neutron in the relevant exit channel
  if(theAngularDistribution!=NULL) 
  {
    theAngularDistribution->SetTarget(theTarget);
    theAngularDistribution->SetNeutron(theNeutron);
  }
  else if(theEnergyAngData!=NULL)
  {
    theEnergyAngData->SetTarget(theTarget);
    theEnergyAngData->SetNeutron(theNeutron);
  }
  
  G4ReactionProductVector * tmpHadrons = NULL;
  G4int i, ii, dummy;
  if(theEnergyAngData != NULL)
  {
    tmpHadrons = theEnergyAngData->Sample(eKinetic);
  }
  else if(theAngularDistribution!= NULL)
  {
    G4bool * Done = new G4bool[nDef];
    for(i=0; i<nDef; i++) Done[i] = false;
    if(tmpHadrons == NULL) 
    {
      tmpHadrons = new G4ReactionProductVector;
    }
    else
    {
      for(i=0; i<tmpHadrons->length(); i++)
      {
	for(ii=0; ii<nDef; ii++)
          if(!Done[ii] && tmpHadrons->at(i)->GetDefinition() == theDefs[ii]) 
              Done[ii] = true;
      }
    }
    G4ReactionProduct * aHadron;
    for(i=0; i<nDef; i++)
    {
      if(!Done[i])
      {
        aHadron = new G4ReactionProduct;
	if(theEnergyDistribution!=NULL)
	{
	  aHadron->SetDefinition(theDefs[i]);
	  aHadron->SetKineticEnergy(theEnergyDistribution->Sample(eKinetic, dummy));
	}
	else
	{
	  G4Exception("No energy distribution to sample from in InelasticBaseFS::BaseApply");
	}
	theAngularDistribution->SampleAndUpdate(*aHadron);
	tmpHadrons->insert(aHadron);
      }
    }
    delete [] Done;
  }
  else
  {
    G4Exception("No data to create the neutrons in NInelasticFS");
  }

  G4ReactionProductVector * thePhotons = NULL;
  if(theFinalStatePhotons!=NULL) 
  {
    // the photon distributions are in the Nucleus rest frame.
    G4ReactionProduct boosted;
    boosted.Lorentz(theNeutron, theTarget);
    G4double anEnergy = boosted.GetKineticEnergy();
    thePhotons = theFinalStatePhotons->GetPhotons(anEnergy);
    for(i=0; i<thePhotons->length(); i++)
    {
      // back to lab
      thePhotons->at(i)->Lorentz(*(thePhotons->at(i)), -1.*theTarget);
    }
  }
  else if(theEnergyAngData!=NULL)
  {
    G4double theGammaEnergy = theEnergyAngData->GetTotalMeanEnergy();
    G4double anEnergy = boosted.GetKineticEnergy();
    theGammaEnergy = anEnergy-theGammaEnergy;
    theGammaEnergy += theNuclearMassDifference;
    G4double eBindProducts = 0;
    G4double eBindN = 0;
    G4double eBindP = 0;
    G4double eBindD = G4NucleiPropertiesTable::GetBindingEnergy(1,2);
    G4double eBindT = G4NucleiPropertiesTable::GetBindingEnergy(1,3);
    G4double eBindHe3 = G4NucleiPropertiesTable::GetBindingEnergy(2,3);
    G4double eBindA = G4NucleiPropertiesTable::GetBindingEnergy(2,4);
    for(i=0; i<tmpHadrons->length(); i++)
    {
      if(tmpHadrons->at(i)->GetDefinition() == G4Neutron::Neutron())
      {
        eBindProducts+=eBindN;
      }
      else if(tmpHadrons->at(i)->GetDefinition() == G4Proton::Proton())
      {
        eBindProducts+=eBindP;
      }
      else if(tmpHadrons->at(i)->GetDefinition() == G4Deuteron::Deuteron())
      {
        eBindProducts+=eBindD;
      }
      else if(tmpHadrons->at(i)->GetDefinition() == G4Triton::Triton())
      {
        eBindProducts+=eBindT;
      }
      else if(tmpHadrons->at(i)->GetDefinition() == G4He3::He3())
      {
        eBindProducts+=eBindHe3;
      }
      else if(tmpHadrons->at(i)->GetDefinition() == G4Alpha::Alpha())
      {
        eBindProducts+=eBindA;
      }
    }
    theGammaEnergy += eBindProducts;
    
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
  }
  
// fill the result
  G4int nSecondaries = tmpHadrons->length();
  G4int nPhotons = 0;
  if(thePhotons!=NULL) nPhotons = thePhotons->length();
  nSecondaries += nPhotons;
  theResult.SetNumberOfSecondaries(nSecondaries);
  G4DynamicParticle * theSec;

  for(i=0; i<nSecondaries-nPhotons; i++)
  {
    theSec = new G4DynamicParticle;    
    theSec->SetDefinition(tmpHadrons->at(i)->GetDefinition());
    theSec->SetMomentum(tmpHadrons->at(i)->GetMomentum());
    theResult.AddSecondary(theSec); 
    delete tmpHadrons->at(i);
  }
  if(thePhotons != NULL)
  {
    for(i=0; i<nPhotons; i++)
    {
      theSec = new G4DynamicParticle;    
      theSec->SetDefinition(thePhotons->at(i)->GetDefinition());
      theSec->SetMomentum(thePhotons->at(i)->GetMomentum());
      theResult.AddSecondary(theSec); 
      delete thePhotons->at(i);
    }
  }
  
// some garbage collection
  delete thePhotons;
  delete tmpHadrons;

// clean up the primary neutron
  theResult.SetStatusChange(fStopAndKill);
}
