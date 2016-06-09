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
//
//
// HPW to implement the choosing of an element for scattering.

#include "G4Types.hh"

#include <fstream>
#include <strstream>
#include <stdlib.h>
#include "G4HadronicProcess.hh"
#include "G4EffectiveCharge.hh"
#include "G4HadProjectile.hh"
#include "G4ElementVector.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Element.hh"
#include "G4ParticleChange.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessManager.hh"
#include "G4StableIsotopes.hh"
#include "G4HadTmpUtil.hh"

#include "G4HadLeadBias.hh"
#include "G4HadronicException.hh"
#include "G4HadReentrentException.hh"

//@@ add model name info, once typeinfo available #include <typeinfo.h>
 
 G4IsoParticleChange * G4HadronicProcess::theIsoResult = NULL;
 G4IsoParticleChange * G4HadronicProcess::theOldIsoResult = NULL;
 G4bool G4HadronicProcess::isoIsEnabled = true;
 
 void G4HadronicProcess::
 EnableIsotopeProductionGlobally()  {isoIsEnabled = true;}
 
 void G4HadronicProcess::
 DisableIsotopeProductionGlobally() {isoIsEnabled = false;}
 
 G4HadronicProcess::G4HadronicProcess( const G4String &processName) :
      G4VDiscreteProcess( processName )
 { 
   isoIsOnAnyway = 0;
   theTotalResult = new G4ParticleChange();
   theCrossSectionDataStore = new G4CrossSectionDataStore();
   aScaleFactor = 1;
   xBiasOn = false;
   if(getenv("SwitchLeadBiasOn")) theBias.push_back(new G4HadLeadBias());
 }

 G4HadronicProcess::~G4HadronicProcess()
 { 
   delete theTotalResult;
   std::for_each(theProductionModels.begin(), 
            theProductionModels.end(), 
	    G4Delete());
   std::for_each(theBias.begin(), 
            theBias.end(), 
	    G4Delete());
 }

 void G4HadronicProcess::RegisterMe( G4HadronicInteraction *a )
 { 
   try{GetManagerPointer()->RegisterMe( a );}   
   catch(G4HadronicException & aE)
   {
     aE.Report(std::cout);
     G4Exception("G4HadronicProcess", "007", FatalException,
                 "Could not register G4HadronicInteraction");
   }
 }

G4double G4HadronicProcess::
GetMeanFreePath(const G4Track &aTrack, G4double, G4ForceCondition *)
{ 
  G4double sigma = 0.0;
  try
  {
    const G4DynamicParticle *aParticle = aTrack.GetDynamicParticle();
    if( !IsApplicable(*aParticle->GetDefinition()))
    {
      G4cout << "Unrecoverable error: "<<G4endl;
      G4ProcessManager * it = aParticle->GetDefinition()->GetProcessManager();
      G4ProcessVector * itv = it->GetProcessList();
      G4cout <<aParticle->GetDefinition()->GetParticleName()<< 
	         " has the following processes:"<<G4endl;
      for(G4int i=0; i<itv->size(); i++)
      {
	G4cout <<"  "<<(*itv)[i]->GetProcessName()<<G4endl;		 
      }
      G4cout << "for kinetic energy "<<aParticle->GetKineticEnergy()<<G4endl;
      G4cout << "and material "<<aTrack.GetMaterial()->GetName()<<G4endl;
      G4Exception("G4HadronicProcess", "007", FatalException,
                   std::string(this->GetProcessName()+
                   " was called for "+
                   aParticle->GetDefinition()->GetParticleName()).c_str() );
    }
    G4Material *aMaterial = aTrack.GetMaterial();
    G4int nElements = aMaterial->GetNumberOfElements();
    
    // returns the mean free path in GEANT4 internal units
    
    const G4double *theAtomicNumDensityVector =
      aMaterial->GetAtomicNumDensityVector();
    
    G4double aTemp = aMaterial->GetTemperature();
        
    for( G4int i=0; i<nElements; ++i )
    {
      G4double xSection =
        GetMicroscopicCrossSection( aParticle, (*aMaterial->GetElementVector())[i], aTemp);
      sigma += theAtomicNumDensityVector[i] * xSection;
    }
    sigma *= aScaleFactor;
    theLastCrossSection = sigma;
  }
  catch(G4HadronicException aR)
  { 
    aR.Report(G4cout); 
    G4Exception("G4HadronicProcess", "007", FatalException,
                "G4HadronicProcess::GetMeanFreePath failed");
  } 
  if( sigma > 0.0 )
    return 1.0/sigma;
  else
    return DBL_MAX;
}

 G4double G4HadronicProcess::GetDistanceToBoundary(const G4Track & aT)
 {
   G4TransportationManager * aTM = 
   G4TransportationManager::GetTransportationManager();
   G4Navigator * aN = aTM->GetNavigatorForTracking();
   G4ThreeVector pGlobalPoint = aT.GetStep()->GetPreStepPoint()->GetPosition();
   G4ThreeVector pDirection = aT.GetMomentumDirection();
   G4double dummy(0);
   G4double result = aN->ComputeStep(pGlobalPoint, pDirection, DBL_MAX, dummy);
   aN->LocateGlobalPointAndSetup(pGlobalPoint);
   return result;
 }

 G4Element * G4HadronicProcess::ChooseAandZ(
  const G4DynamicParticle *aParticle, const G4Material *aMaterial )
  {
    static G4bool noIsotopeWiseCrossSections=getenv("GHAD_DISABLE_ISOTOPE_WISE_CROSS_SECTIONS");
    static G4StableIsotopes theIso;
    currentZ = 0;
    currentN = 0;
    const G4int numberOfElements = aMaterial->GetNumberOfElements();
    const G4ElementVector *theElementVector = aMaterial->GetElementVector();
    G4int i;
    if( numberOfElements == 1 ) 
    {
      currentZ = G4double( ((*theElementVector)[0])->GetZ());
      G4int localZ = G4lrint(currentZ);
      if(noIsotopeWiseCrossSections)
      {
        currentN = (*theElementVector)[0]->GetN();
      }
      else
      {
        G4double * running = new G4double[theIso.GetNumberOfIsotopes(localZ)];
	for (i=0; i<theIso.GetNumberOfIsotopes(localZ); i++)
        {
          G4double fracInPercent=theIso.GetAbundance(theIso.GetFirstIsotope(localZ)+i);
          G4double runningA=theIso.GetIsotopeNucleonCount(theIso.GetFirstIsotope(localZ)+i);
	  running[i]=fracInPercent*pow(runningA, 2./3.);
	  // rough approximation; to get it better, redesign getMSC to not use G4Element, see also below
	  if(i!=0) running[i] += running[i-1];
        }
	G4double trial = G4UniformRand();
	G4double sum = running[theIso.GetNumberOfIsotopes(localZ)-1];
	for(i=0; i<theIso.GetNumberOfIsotopes(localZ); i++)
	{
	  currentN = theIso.GetIsotopeNucleonCount(theIso.GetFirstIsotope(localZ)+i);
	  if(running[i]/sum>trial) break;
	}
	delete [] running;
      }
      targetNucleus.SetParameters(currentN, currentZ);
      return (*theElementVector)[0];
    }
    
    const G4double *theAtomicNumberDensity = aMaterial->GetAtomicNumDensityVector();
    G4double aTemp = aMaterial->GetTemperature();
    G4double crossSectionTotal = 0;
    std::vector<G4double> runningSum;
    for( i=0; i < numberOfElements; ++i )
    {
      runningSum.push_back(theAtomicNumberDensity[i] *
        dispatch->GetMicroscopicCrossSection( aParticle, (*theElementVector)[i], aTemp));
      crossSectionTotal+=runningSum[i];
    }
    
    G4double random = G4UniformRand();
    for( i=0; i < numberOfElements; ++i )
    { 
      if(i!=0) runningSum[i]+=runningSum[i-1];
      if( random<=runningSum[i]/crossSectionTotal )
      {
        currentZ = G4double( ((*theElementVector)[i])->GetZ());
        G4int localZ = G4lrint(currentZ);
        if(noIsotopeWiseCrossSections)
        {
          currentN = ((*theElementVector)[i])->GetN();
        }
	else
	{
          G4double * running = new G4double[theIso.GetNumberOfIsotopes(localZ)];
	  for (i=0; i<theIso.GetNumberOfIsotopes(localZ); i++)
          {
            G4double fracInPercent=theIso.GetAbundance(theIso.GetFirstIsotope(localZ)+i);
            G4double runningA=theIso.GetIsotopeNucleonCount(theIso.GetFirstIsotope(localZ)+i);
	    running[i]=fracInPercent*pow(runningA, 2./3.);
	    if(i!=0) running[i] += running[i-1];
          }
	  G4double trial = G4UniformRand();
	  for(i=0; i<theIso.GetNumberOfIsotopes(localZ); i++)
	  {
	    currentN = theIso.GetIsotopeNucleonCount(theIso.GetFirstIsotope(localZ)+i);
	    if(running[i]/running[theIso.GetNumberOfIsotopes(localZ)-1]>trial) break;
	  }
  	  delete [] running;
	}
	targetNucleus.SetParameters(currentN, currentZ);
        return (*theElementVector)[i];
      }
    }
    currentZ = G4double((*theElementVector)[numberOfElements-1]->GetZ());
    G4int localZ = G4lrint(currentZ);
    if(noIsotopeWiseCrossSections)
    {
      currentN = (*theElementVector)[numberOfElements-1]->GetN();
    }
    else
    {
      G4double * running = new G4double[theIso.GetNumberOfIsotopes(localZ)];
      for (i=0; i<theIso.GetNumberOfIsotopes(localZ); i++)
      {
        G4double fracInPercent=theIso.GetAbundance(theIso.GetFirstIsotope(localZ)+i);
        G4double runningA=theIso.GetIsotopeNucleonCount(theIso.GetFirstIsotope(localZ)+i);
	running[i]=fracInPercent*pow(runningA, 2./3.);
	// rough approximation; to get it better, redesign getMSC to not use G4Element
	if(i!=0) running[i] += running[i-1];
      }
      G4double trial = G4UniformRand();
      for(i=0; i<theIso.GetNumberOfIsotopes(localZ); i++)
      {
	currentN = theIso.GetIsotopeNucleonCount(theIso.GetFirstIsotope(localZ)+i);
	if(running[i]/running[theIso.GetNumberOfIsotopes(localZ)-1]>trial) break;
      }
      delete [] running;
    }
    targetNucleus.SetParameters(currentN, currentZ);
    return (*theElementVector)[numberOfElements-1];
  }
 
 G4VParticleChange *G4HadronicProcess::GeneralPostStepDoIt(
  const G4Track &aTrack, const G4Step &)
  {
    // G4cout << theNumberOfInteractionLengthLeft<<G4endl;
    const G4DynamicParticle *aParticle = aTrack.GetDynamicParticle();
    G4Material *aMaterial = aTrack.GetMaterial();
    G4double originalEnergy = aParticle->GetKineticEnergy();
    G4double kineticEnergy = originalEnergy;
    if(aParticle->GetDefinition()->GetBaryonNumber()>1.5)
    {
      kineticEnergy/=aParticle->GetDefinition()->GetBaryonNumber();
    }
    G4Element * anElement = 0;
    try
    {
      anElement = ChooseAandZ( aParticle, aMaterial );
    }
    catch(G4HadronicException & aR)
    {
      aR.Report(G4cout);
      G4cout << "Unrecoverable error for:"<<G4endl;
      G4cout << " - Particle energy[GeV] = "<< originalEnergy/GeV<<G4endl;
      G4cout << " - Material = "<<aMaterial->GetName()<<G4endl;
      G4cout << " - Particle type = "
             <<aParticle->GetDefinition()->GetParticleName()<<G4endl;
      G4Exception("G4HadronicProcess", "007", FatalException,
                  "GeneralPostStepDoIt failed on element selection.");
    }
    try
    {
    theInteraction = ChooseHadronicInteraction( kineticEnergy,
                                                aMaterial, anElement );
    }
    catch(G4HadronicException & aE)
    {
      aE.Report(std::cout);
      G4cout << "Unrecoverable error for:"<<G4endl;
      G4cout << " - Particle energy[GeV] = "<< originalEnergy/GeV<<G4endl;
      G4cout << " - Material = "<<aMaterial->GetName()<<G4endl;
      G4cout << " - Particle type = "<<aParticle->GetDefinition()->GetParticleName()<<G4endl;
      G4Exception("G4HadronicProcess", "007", FatalException,
                  "ChooseHadronicInteraction failed.");
    }
    G4HadProjectile thePro(aTrack);
    
    G4HadFinalState *result = 0;
    G4int reentryCount = 0;
    do
    {
      try
      {
         result = theInteraction->ApplyYourself( thePro, targetNucleus);
      }
      catch(G4HadReentrentException aR)
      {
        aR.Report(G4cout);
	G4cout << " G4HadronicProcess re-entering the ApplyYourself call for"<<G4endl;
        G4cout << " - Particle energy[GeV] = "<< originalEnergy/GeV<<G4endl;
        G4cout << " - Material = "<<aMaterial->GetName()<<G4endl;
        G4cout << " - Particle type = "<<aParticle->GetDefinition()->GetParticleName()<<G4endl;
	result = 0; // here would still be leaking...
	if(reentryCount>100)
	{
           G4Exception("G4HadronicProcess", "007", FatalException,
	  "GetHadronicProcess: Reentering ApplyYourself too often - GeneralPostStepDoIt failed.");  
	}
        G4Exception("G4HadronicProcess", "007", FatalException,
	  "GetHadronicProcess: GeneralPostStepDoIt failed (Reentering ApplyYourself not yet supported.)");  
      }
      catch(G4HadronicException aR)
      {
        aR.Report(G4cout);
	G4cout << " G4HadronicProcess failed in ApplyYourself call for"<<G4endl;
        G4cout << " - Particle energy[GeV] = "<< originalEnergy/GeV<<G4endl;
        G4cout << " - Material = "<<aMaterial->GetName()<<G4endl;
        G4cout << " - Particle type = "<<aParticle->GetDefinition()->GetParticleName()<<G4endl;
        G4Exception("G4HadronicProcess", "007", FatalException,
	            "GeneralPostStepDoIt failed.");
      }
    }
    while(!result);
    if(result->GetStatusChange() == isAlive && thePro.GetDefinition() != aTrack.GetDefinition())
    {
      G4DynamicParticle * aP = const_cast<G4DynamicParticle *>(aTrack.GetDynamicParticle());
      aP->SetDefinition(const_cast<G4ParticleDefinition *>(thePro.GetDefinition()));
    }
    result->SetTrafoToLab(thePro.GetTrafoToLab());
    for(G4int i=0; i<result->GetNumberOfSecondaries(); i++)
    {
      G4DynamicParticle* aSecTrack = result->GetSecondary(i)->GetParticle();
      if(aSecTrack->GetDefinition()->GetPDGCharge()>1.5)
      {
         G4EffectiveCharge aCalculator;
	 G4double charge = aCalculator.GetCharge(aMaterial, aSecTrack->GetKineticEnergy(),
	                                        aSecTrack->GetDefinition()->GetPDGMass(),
						aSecTrack->GetDefinition()->GetPDGCharge());
         if(getenv("GHADChargeDebug")) 
         {
            std::cout << "Recoil fractional charge is "
                      << charge/aSecTrack->GetDefinition()->GetPDGCharge()<<" "
                      << charge <<" "<<aSecTrack->GetDefinition()->GetPDGCharge()<<std::endl;
         }
	 aSecTrack->SetCharge(charge);
      }
    }

    if(getenv("HadronicDoitLogging") )
    {
      G4cout << "HadronicDoitLogging "
             << GetProcessName() <<" "
             << aParticle->GetDefinition()->GetPDGEncoding()<<" "
	     << originalEnergy<<" "
	     << aParticle->GetMomentum()<<" "
	     << targetNucleus.GetN()<<" "
	     << targetNucleus.GetZ()<<" "
	     << G4endl;
    }
    ClearNumberOfInteractionLengthLeft();
    if(isoIsOnAnyway!=-1)
    {
      if(isoIsEnabled||isoIsOnAnyway)
      {
        result = DoIsotopeCounting(result, aTrack, targetNucleus);
      }
    }
    
    G4double e=aTrack.GetKineticEnergy();
    if(e<5*GeV)
    {
      for(size_t i=0; i<theBias.size(); i++)
      {
        result = theBias[i]->Bias(result);
      }
    }
    FillTotalResult(result, aTrack);
    return theTotalResult;
  }

  G4HadFinalState * G4HadronicProcess::
  DoIsotopeCounting(G4HadFinalState * aResult,
                    const G4Track & aTrack,
                    const G4Nucleus & aNucleus)
  {
    // get the PC from iso-production
    if(theOldIsoResult) delete theOldIsoResult;
    if(theIsoResult) delete theIsoResult;
    theIsoResult = new G4IsoParticleChange;
    G4bool done = false;
    G4IsoResult * anIsoResult = NULL;
    for(unsigned int i=0; i<theProductionModels.size(); i++)
    {
      anIsoResult = theProductionModels[i]->GetIsotope(aTrack, aNucleus);
      if(anIsoResult!=NULL)
      {
        done = true;
        break;
      }
    }
    // if none in charge, use default iso production
    if(!done) anIsoResult = ExtractResidualNucleus(aTrack, aNucleus, aResult); 
        
    // Add all info explicitely and add typename from model called.
    theIsoResult->SetIsotope(anIsoResult->GetIsotope());
    theIsoResult->SetProductionPosition(aTrack.GetPosition());
    theIsoResult->SetProductionTime(aTrack.GetGlobalTime());
    theIsoResult->SetParentParticle(*aTrack.GetDynamicParticle());
    theIsoResult->SetMotherNucleus(anIsoResult->GetMotherNucleus());
//    theIsoResult->SetProducer(typeid(*theInteraction).name()); @@@@@@@
    G4String aWorkaround("WaitingForTypeidToBeAvailableInCompilers"); // @@@@@  workaround for DEC.
    theIsoResult->SetProducer(aWorkaround);
    
    delete anIsoResult;

    return aResult;
  }
  
  G4IsoResult * G4HadronicProcess::
  ExtractResidualNucleus(const G4Track & ,
                         const G4Nucleus & aNucleus,
                         G4HadFinalState * aResult)
  {
    G4double A = aNucleus.GetN();
    G4double Z = aNucleus.GetZ();
    G4double bufferA = 0;
    G4double bufferZ = 0;
    
    // loop over aResult, and decrement A, Z accordingly
    // cash the max
    for(G4int i=0; i<aResult->GetNumberOfSecondaries(); i++)
    {
      G4HadSecondary* aSecTrack = aResult->GetSecondary(i);
      if(bufferA<aSecTrack->GetParticle()->GetDefinition()->GetBaryonNumber())
      {
        bufferA = aSecTrack->GetParticle()->GetDefinition()->GetBaryonNumber();
        bufferZ = aSecTrack->GetParticle()->GetDefinition()->GetPDGCharge();
      }
      Z-=aSecTrack->GetParticle()->GetDefinition()->GetPDGCharge();
      A-=aSecTrack->GetParticle()->GetDefinition()->GetBaryonNumber();
    }
    
    // if the fragment was part of the final state, it is 
    // assumed to be the heaviest secondary.
    if(A<0.1)
    {
       A = bufferA;
       Z = bufferZ;
    }
    
    // prepare the IsoResult.
    char the1[100] = {""};
    std::ostrstream ost1(the1, 100, std::ios::out);
    ost1 <<Z<<"_"<<A<<"\0";
    G4String * biff = new G4String(the1);
    G4IsoResult * theResult = new G4IsoResult(*biff, aNucleus);
    
    // cleaning up.
    delete biff;
    
    return theResult;
  }

G4double G4HadronicProcess::
XBiasSurvivalProbability()
{
  G4double result = 0;
  G4double nLTraversed = GetTotalNumberOfInteractionLengthTraversed();
  G4double biasedProbability = 1.-exp(-nLTraversed);
  G4double realProbability = 1-exp(-nLTraversed/aScaleFactor);
  result = (biasedProbability-realProbability)/biasedProbability;
  return result;
}

G4double G4HadronicProcess::
XBiasSecondaryWeight()
{
  G4double result = 0;
  G4double nLTraversed = GetTotalNumberOfInteractionLengthTraversed();
  result = 1./aScaleFactor*exp(-nLTraversed/aScaleFactor*(1-1./aScaleFactor));
  return result;
}

void G4HadronicProcess::FillTotalResult(G4HadFinalState * aR, const G4Track & aT)
{
//          G4cout << "############# Entry debug "
//	         <<GetProcessName()<<" "
//		 <<aT.GetDynamicParticle()->GetDefinition()->GetParticleName()<<" "
//		 <<aT.GetDynamicParticle()<<" "
 //                <<aScaleFactor<<" "
//		 <<aT.GetWeight()<<" "
//		 <<G4endl;
      theTotalResult->Clear();
      theTotalResult->SetLocalEnergyDeposit(0.);
      theTotalResult->Initialize(aT);
      theTotalResult->SetSecondaryWeightByProcess(true);
      theTotalResult->SetStatusChange(fAlive);
      G4double rotation = 2.*pi*G4UniformRand();
      G4ThreeVector it(0., 0., 1.);
      /*
      if(xBiasOn)
      {
        G4cout << "BiasDebug "<<GetProcessName()<<" "
                              <<aScaleFactor<<" "
                              <<XBiasSurvivalProbability()<<" "
		  	      <<XBiasSecondaryWeight()<<" "
			      <<G4endl;
      }
      */
      if(aR->GetStatusChange()==stopAndKill)
      {
	if( xBiasOn && G4UniformRand()<XBiasSurvivalProbability() )
	{
 	  theTotalResult->SetWeightChange( XBiasSurvivalProbability()*aT.GetWeight() );
	}
	else
	{
          theTotalResult->SetStatusChange(fStopAndKill);
          theTotalResult->SetEnergyChange( 0.0 );
	}
      }
      else if(aR->GetStatusChange()!=stopAndKill )
      {
        if(aR->GetStatusChange()==suspend)
        {
	  if(xBiasOn)
	  {
            G4Exception("G4HadronicProcess", "007", FatalException,
	                "Cannot cross-section bias a process that suspends tracks.");
	  }
          theTotalResult->SetStatusChange(fSuspend);
        }
	if(xBiasOn && G4UniformRand()<XBiasSurvivalProbability())
	{
 	  theTotalResult->SetWeightChange( XBiasSurvivalProbability()*aT.GetWeight() );
	  G4double newWeight = aR->GetWeightChange()*aT.GetWeight();
          G4DynamicParticle * aNew = new G4DynamicParticle(aT.GetDefinition(),
	                                                   aR->GetEnergyChange(),
							   aR->GetMomentumChange());
	  G4HadSecondary * theSec = new G4HadSecondary(aNew, newWeight);
	  aR->AddSecondary(theSec);
	}
	else
	{
	  G4double newWeight = aR->GetWeightChange()*aT.GetWeight();
	  theTotalResult->SetWeightChange(newWeight); // This is multiplicative
	  if(aR->GetEnergyChange()>-.5) theTotalResult->SetEnergyChange(aR->GetEnergyChange());
	  G4LorentzVector newDirection(aR->GetMomentumChange().unit(), 1.);
	  newDirection*=aR->GetTrafoToLab();
	  theTotalResult->SetMomentumDirectionChange(newDirection.vect());
	}
      }
      else
      {
          G4cerr << "Track status is "<< aR->GetStatusChange()<<G4endl;
	  G4Exception("G4HadronicProcess", "007", FatalException,
	              "use of unsupported track-status.");
      }

      theTotalResult->SetLocalEnergyDeposit(aR->GetLocalEnergyDeposit());
      theTotalResult->SetNumberOfSecondaries(aR->GetNumberOfSecondaries());

      for(G4int i=0; i<aR->GetNumberOfSecondaries(); i++)
      {
        G4LorentzVector theM = aR->GetSecondary(i)->GetParticle()->Get4Momentum();
        theM.rotate(rotation, it);
	theM*=aR->GetTrafoToLab();
	aR->GetSecondary(i)->GetParticle()->Set4Momentum(theM);
	G4double time = aR->GetSecondary(i)->GetTime();
	if(time<0) time = aT.GetGlobalTime();
        G4Track* track = new G4Track(aR->GetSecondary(i)->GetParticle(),
				     aT.GetGlobalTime(),
				     aT.GetPosition());
	G4double newWeight = aT.GetWeight()*aR->GetSecondary(i)->GetWeight();
	//static G4double pinelcount=0;
	if(xBiasOn) newWeight *= XBiasSecondaryWeight();
        /*  G4cout << "#### ParticleDebug "
	         <<GetProcessName()<<" "
		 <<aR->GetSecondary(i)->GetParticle()->GetDefinition()->GetParticleName()<<" "
                 <<aScaleFactor<<" "
                 <<XBiasSurvivalProbability()<<" "
	  	 <<XBiasSecondaryWeight()<<" "
		 <<aT.GetWeight()<<" "
		 <<aR->GetSecondary(i)->GetWeight()<<" "
		 <<aR->GetSecondary(i)->GetParticle()<<" "
		 <<G4endl;*/
	track->SetWeight(newWeight);
	G4double trackDeb = track->GetKineticEnergy();
	if( (  trackDeb<0 
	   || (trackDeb>aT.GetKineticEnergy()+1*GeV) ) && getenv("GHADEnergyBalanceDebug") )
	{
	  G4cout << "Debugging hadronic processes: "<<track->GetKineticEnergy()
	         <<" "<<aT.GetKineticEnergy()
		 <<" "<<GetProcessName()
		 <<" "<<aT.GetDefinition()->GetParticleName() 
		 <<G4endl;
	}
	/*if(GetProcessName()=="PhotonInelastic")
	{
	  if(aR->GetSecondary(i)->GetParticle()->GetDefinition()==G4Neutron::NeutronDefinition())
	  {
	   pinelcount+= newWeight;
	   G4cout << "=======> Neutrons from gamma-nuclear "<<pinelcount<<G4endl;
	  }
	}*/
	
	theTotalResult->AddSecondary(track);
      }
      aR->Clear();
      return;
}
 /* end of file */
