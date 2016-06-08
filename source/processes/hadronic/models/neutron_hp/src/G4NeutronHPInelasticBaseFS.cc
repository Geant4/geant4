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
#include "G4NeutronHPInelasticBaseFS.hh"
#include "G4Nucleus.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4Electron.hh"
#include "G4NeutronHPDataUsed.hh"

void G4NeutronHPInelasticBaseFS::InitGammas(G4double AR, G4double ZR)
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
    
   G4double eps = 0.001;
   theNuclearMassDifference = 
       G4NucleiPropertiesTable::GetBindingEnergy(static_cast<G4int>(ZR+eps),static_cast<G4int>(AR+eps)) -
       G4NucleiPropertiesTable::GetBindingEnergy(static_cast<G4int>(theBaseZ+eps), static_cast<G4int>(theBaseA+eps));
   theGammas.Init(theGammaData);
   delete aName;
}

void G4NeutronHPInelasticBaseFS::Init (G4double A, G4double Z, G4String & dirName, G4String & bit)
{
  gammaPath = "/Inelastic/Gammas/";
    if(!getenv("NeutronHPCrossSections")) 
       G4Exception("Please setenv NeutronHPCrossSections to point to the neutron cross-section files.");
  G4String tBase = getenv("NeutronHPCrossSections");
  gammaPath = tBase+gammaPath;
  G4String tString = dirName;
  G4bool dbool;
  G4NeutronHPDataUsed aFile = theNames.GetName(static_cast<G4int>(A), static_cast<G4int>(Z), tString, bit, dbool);
  G4String filename = aFile.GetName();
  theBaseA = aFile.GetA();
  theBaseZ = aFile.GetZ();
  if(!dbool || ( Z<2.5 && ( abs(theBaseZ - Z)>0.0001 || abs(theBaseA - A)>0.0001)))
  {
    if(getenv("NeutronHPNamesLogging")) G4cout << "Skipped = "<< filename <<" "<<A<<" "<<Z<<G4endl;
    hasAnyData = false;
    hasFSData = false; 
    hasXsec = false;
    return;
  }
    theBaseA = A;
    theBaseZ = G4int(Z+.5);
#ifdef G4USE_STD_NAMESPACE
  G4std::ifstream theData(filename, G4std::ios::in);
#else
  G4std::ifstream theData(filename, G4std::ios::in|G4std::ios::nocreate);
#endif
  if(!(theData))
  {
    hasAnyData = false;
    hasFSData = false; 
    hasXsec = false;
    theData.close();
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
  theData.close();
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
  targetMass = ( G4NucleiPropertiesTable::GetAtomicMass(static_cast<G4int>(theBaseZ+eps), static_cast<G4int>(theBaseA+eps))-
                            theBaseZ*G4Electron::ElectronDefinition()->GetPDGMass() ) /
               G4Neutron::Neutron()->GetPDGMass();
  if(theEnergyAngData!=NULL)
      targetMass = theEnergyAngData->GetTargetMass();
  if(theAngularDistribution!=NULL)
      targetMass = theAngularDistribution->GetTargetMass();
  G4Nucleus aNucleus;
  G4ReactionProduct theTarget; 
  G4ThreeVector neuVelo = (1./incidentParticle->GetDefinition()->GetPDGMass())*theNeutron.GetMomentum();
  theTarget = aNucleus.GetBiasedThermalNucleus( targetMass, neuVelo, theTrack.GetMaterial()->GetTemperature());

// prepare energy in target rest frame
  G4ReactionProduct boosted;
  boosted.Lorentz(theNeutron, theTarget);
  eKinetic = boosted.GetKineticEnergy();
  G4double orgMomentum = boosted.GetMomentum().mag();
  
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
  G4int ii, dummy;
  unsigned int i;
  if(theEnergyAngData != NULL)
  {
    tmpHadrons = theEnergyAngData->Sample(eKinetic);
  }
  else if(theAngularDistribution!= NULL)
  {
    G4bool * Done = new G4bool[nDef];
    G4int i0;
    for(i0=0; i0<nDef; i0++) Done[i0] = false;
    if(tmpHadrons == NULL) 
    {
      tmpHadrons = new G4ReactionProductVector;
    }
    else
    {
      for(i=0; i<tmpHadrons->size(); i++)
      {
	for(ii=0; ii<nDef; ii++)
          if(!Done[ii] && tmpHadrons->operator[](i)->GetDefinition() == theDefs[ii]) 
              Done[ii] = true;
      }
    }
    G4ReactionProduct * aHadron;
    G4double localMass = ( G4NucleiPropertiesTable::GetAtomicMass(static_cast<G4int>(theBaseZ+eps), static_cast<G4int>(theBaseA+eps))-
                            theBaseZ*G4Electron::ElectronDefinition()->GetPDGMass() );
    G4ThreeVector bufferedDirection(0,0,0);
    for(i0=0; i0<nDef; i0++)
    {
      if(!Done[i0])
      {
        aHadron = new G4ReactionProduct;
	if(theEnergyDistribution!=NULL)
	{
	  aHadron->SetDefinition(theDefs[i0]);
	  aHadron->SetKineticEnergy(theEnergyDistribution->Sample(eKinetic, dummy));
	}
	else if(nDef == 1)
	{
	  aHadron->SetDefinition(theDefs[i0]);
	  aHadron->SetKineticEnergy(eKinetic);
	}
	else if(nDef == 2)
	{
	  aHadron->SetDefinition(theDefs[i0]);
	  aHadron->SetKineticEnergy(50*MeV);
	}
	else
	{
	  G4Exception("No energy distribution to sample from in InelasticBaseFS::BaseApply");
	}
	theAngularDistribution->SampleAndUpdate(*aHadron);
	if(theEnergyDistribution==NULL && nDef == 2)
	{
	  if(i0==0)
	  {
	    G4double m1 = theDefs[0]->GetPDGMass();
	    G4double m2 = theDefs[1]->GetPDGMass();
	    G4double mn = G4Neutron::Neutron()->GetPDGMass();
	    G4int z1 = static_cast<G4int>(theBaseZ+eps-theDefs[0]->GetPDGCharge()-theDefs[1]->GetPDGCharge());
	    G4int a1 = static_cast<G4int>(theBaseA+eps)-theDefs[0]->GetBaryonNumber()-theDefs[1]->GetBaryonNumber();
	    G4double concreteMass = G4NucleiPropertiesTable::GetAtomicMass(z1, a1)-z1*G4Electron::ElectronDefinition()->GetPDGMass();
	    G4double availableEnergy = eKinetic+mn+localMass-m1-m2-concreteMass;
	    // available kinetic energy in CMS (non relativistic)
	    G4double emin = availableEnergy+m1+m2 - sqrt((m1+m2)*(m1+m2)+orgMomentum*orgMomentum);
	    G4double p1=sqrt(2.*m2*emin);
	    bufferedDirection = p1*aHadron->GetMomentum().unit();
	    if(getenv("HTOKEN")) // @@@@@ verify the nucleon counting...
	    { 
	      G4cout << "HTOKEN "<<z1<<" "<<theBaseZ<<" "<<a1<<" "<<theBaseA<<" "<<availableEnergy<<" "
	             << emin<<G4endl;
            }
	  }
	  else
	  {
	    bufferedDirection = -bufferedDirection;
	  }
	  // boost from cms to lab
	  if(getenv("HTOKEN")) 
	  {
	    G4cout << " HTOKEN "<<bufferedDirection.mag2()<<G4endl;
	  }
	  aHadron->SetTotalEnergy( sqrt(aHadron->GetMass()*aHadron->GetMass()
	                              +bufferedDirection.mag2()) );
	  aHadron->SetMomentum(bufferedDirection);
          aHadron->Lorentz(*aHadron, -1.*(theTarget+theNeutron)); 
	  if(getenv("HTOKEN")) 
	  {
	    G4cout << "  HTOKEN "<<aHadron->GetTotalEnergy()<<" "<<aHadron->GetMomentum()<<G4endl;
	  }
	}
	tmpHadrons->push_back(aHadron);
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
    if(thePhotons!=NULL)
    {
      for(i=0; i<thePhotons->size(); i++)
      {
        // back to lab
        thePhotons->operator[](i)->Lorentz(*(thePhotons->operator[](i)), -1.*theTarget);
      }
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
    for(i=0; i<tmpHadrons->size(); i++)
    {
      if(tmpHadrons->operator[](i)->GetDefinition() == G4Neutron::Neutron())
      {
        eBindProducts+=eBindN;
      }
      else if(tmpHadrons->operator[](i)->GetDefinition() == G4Proton::Proton())
      {
        eBindProducts+=eBindP;
      }
      else if(tmpHadrons->operator[](i)->GetDefinition() == G4Deuteron::Deuteron())
      {
        eBindProducts+=eBindD;
      }
      else if(tmpHadrons->operator[](i)->GetDefinition() == G4Triton::Triton())
      {
        eBindProducts+=eBindT;
      }
      else if(tmpHadrons->operator[](i)->GetDefinition() == G4He3::He3())
      {
        eBindProducts+=eBindHe3;
      }
      else if(tmpHadrons->operator[](i)->GetDefinition() == G4Alpha::Alpha())
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
        for(unsigned int ii=0; ii<theOtherPhotons->size(); ii++)
        {
          thePhotons->push_back(theOtherPhotons->operator[](ii));
        }
        delete theOtherPhotons; 
      }
      theGammaEnergy -= theGammas.GetLevelEnergy(iLevel);
      if(iLevel == -1) break;
    }
  }
  
// fill the result
  unsigned int nSecondaries = tmpHadrons->size();
  unsigned int nPhotons = 0;
  if(thePhotons!=NULL) nPhotons = thePhotons->size();
  nSecondaries += nPhotons;
  theResult.SetNumberOfSecondaries(nSecondaries);
  G4DynamicParticle * theSec;

  for(i=0; i<nSecondaries-nPhotons; i++)
  {
    theSec = new G4DynamicParticle;    
    theSec->SetDefinition(tmpHadrons->operator[](i)->GetDefinition());
    theSec->SetMomentum(tmpHadrons->operator[](i)->GetMomentum());
    theResult.AddSecondary(theSec); 
    delete tmpHadrons->operator[](i);
  }
  if(thePhotons != NULL)
  {
    for(i=0; i<nPhotons; i++)
    {
      theSec = new G4DynamicParticle;    
      theSec->SetDefinition(thePhotons->operator[](i)->GetDefinition());
      theSec->SetMomentum(thePhotons->operator[](i)->GetMomentum());
      theResult.AddSecondary(theSec); 
      delete thePhotons->operator[](i);
    }
  }
  
// some garbage collection
  delete thePhotons;
  delete tmpHadrons;

// clean up the primary neutron
  theResult.SetStatusChange(fStopAndKill);
}
