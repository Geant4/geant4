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
#include "G4NoModelFound.hh"
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
 }

 G4HadronicProcess::~G4HadronicProcess()
 { 
   delete theTotalResult;
   for_each(theProductionModels.begin(), 
            theProductionModels.end(), 
	    Delete<G4VIsotopeProduction>());
   for_each(theBias.begin(), 
            theBias.end(), 
	    Delete<G4VLeadingParticleBiasing>());
 }

 void G4HadronicProcess::RegisterMe( G4HadronicInteraction *a )
 { GetManagerPointer()->RegisterMe( a ); }

 G4double G4HadronicProcess::
 GetMeanFreePath(const G4Track &aTrack, G4double, G4ForceCondition *)
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
      G4Exception( this->GetProcessName()+
                   " was called for "+
                   aParticle->GetDefinition()->GetParticleName() );
    }
    G4Material *aMaterial = aTrack.GetMaterial();
    G4int nElements = aMaterial->GetNumberOfElements();
    
    // returns the mean free path in GEANT4 internal units
    
    const G4double *theAtomicNumDensityVector =
      aMaterial->GetAtomicNumDensityVector();
    
    G4double aTemp = aMaterial->GetTemperature();
        
    G4double sigma = 0.0;
    for( G4int i=0; i<nElements; ++i )
    {
      G4double xSection =
        GetMicroscopicCrossSection( aParticle, (*aMaterial->GetElementVector())[i], aTemp);
      sigma += theAtomicNumDensityVector[i] * xSection;
    }
    sigma *= aScaleFactor;
    theLastCrossSection = sigma;
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
    currentZ = 0;
    currentN = 0;
    const G4int numberOfElements = aMaterial->GetNumberOfElements();
    const G4ElementVector *theElementVector = aMaterial->GetElementVector();
    
    if( numberOfElements == 1 ) 
    {
      currentZ = G4double( ((*theElementVector)[0])->GetZ());
      currentN = (*theElementVector)[0]->GetN();
      targetNucleus.SetParameters(currentN, currentZ);
      return (*theElementVector)[0];
    }
    
    const G4double *theAtomicNumberDensity = aMaterial->GetAtomicNumDensityVector();
    G4double aTemp = aMaterial->GetTemperature();
    G4double crossSectionTotal = 0;
    G4int i;
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
      if( random<=runningSum[i]/crossSectionTotal )
      {
        currentZ = G4double( ((*theElementVector)[i])->GetZ());
        currentN = ((*theElementVector)[i])->GetN();
        targetNucleus.SetParameters(currentN, currentZ);
        return (*theElementVector)[i];
      }
    }
    currentZ = G4double((*theElementVector)[numberOfElements-1]->GetZ());
    currentN = (*theElementVector)[numberOfElements-1]->GetN();
    targetNucleus.SetParameters(currentN, currentZ);
    return (*theElementVector)[numberOfElements-1];
  }
 
 G4VParticleChange *G4HadronicProcess::GeneralPostStepDoIt(
  const G4Track &aTrack, const G4Step &)
  {
    // G4cout << theNumberOfInteractionLengthLeft<<G4endl;
    const G4DynamicParticle *aParticle = aTrack.GetDynamicParticle();
    G4Material *aMaterial = aTrack.GetMaterial();
    G4double kineticEnergy = aParticle->GetKineticEnergy();
    G4Element * anElement = ChooseAandZ( aParticle, aMaterial );
    try
    {
    theInteraction = ChooseHadronicInteraction( kineticEnergy,
                                                aMaterial, anElement );
    }
    catch(G4NoModelFound * it)
    {
      delete it;
      G4cout << "Unrecoverable error for:"<<G4endl;
      G4cout << " - Particle energy[GeV] = "<< kineticEnergy/GeV<<G4endl;
      G4cout << " - Material = "<<aMaterial->GetName()<<G4endl;
      G4cout << " - Particle type = "<<aParticle->GetDefinition()->GetParticleName()<<G4endl;
      G4Exception("GetHadronicProcess: No model found for this energy range");
    }
    G4HadProjectile thePro(aTrack);
    G4HadFinalState *result =
      theInteraction->ApplyYourself( thePro, targetNucleus);
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
	 G4double charge = aCalculator.GetCharge(aMaterial, kineticEnergy,
	                                        aSecTrack->GetDefinition()->GetPDGMass(),
						aSecTrack->GetDefinition()->GetPDGCharge());
	 aSecTrack->SetCharge(charge);
      }
    }

    if(getenv("HadronicDoitLogging") )
    {
      G4cout << "HadronicDoitLogging "
             << GetProcessName() <<" "
             << aParticle->GetDefinition()->GetPDGEncoding()<<" "
	     << kineticEnergy<<" "
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
    
    for(size_t i=0; i<theBias.size(); i++)
    {
      result = theBias[i]->Bias(result);
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
  G4double biasedProbability = exp(-nLTraversed);
  G4double realProbability = 1./aScaleFactor*exp(-nLTraversed/aScaleFactor);
  result = (biasedProbability-realProbability)/biasedProbability;
  return result;
}

G4double G4HadronicProcess::
XBiasSecondaryWeight()
{
  G4double result = 0;
  G4double nLTraversed = GetTotalNumberOfInteractionLengthTraversed();
  G4double biasedProbability = exp(-nLTraversed);
  G4double realProbability = 1./aScaleFactor*exp(-nLTraversed/aScaleFactor);
  result = realProbability/biasedProbability;
  return result;
}

void G4HadronicProcess::FillTotalResult(G4HadFinalState * aR, const G4Track & aT)
{
      theTotalResult->Clear();
      theTotalResult->Initialize(aT);
      theTotalResult->SetSecondaryWeightByProcess(true);
      theTotalResult->SetStatusChange(fAlive);
      G4double rotation = 2.*pi*G4UniformRand();
      G4ThreeVector it(0., 0., 1.);
      if(xBiasOn)
      {
        G4cout << "BiasDebug "<<GetProcessName()<<" "
                              <<aScaleFactor<<" "
                              <<XBiasSurvivalProbability()<<" "
		  	      <<XBiasSecondaryWeight()<<" "
			      <<G4endl;
      }
      if(aR->GetStatusChange()==stopAndKill)
      {
	if( xBiasOn && G4UniformRand()<XBiasSurvivalProbability() )
	{
 	  theTotalResult->SetWeightChange( XBiasSurvivalProbability() );
	}
	else
	{
          theTotalResult->SetStatusChange(fStopAndKill);
          theTotalResult->SetEnergyChange( 0.0 );
	}
      }
      if(aR->GetStatusChange()==suspend)
      {
        theTotalResult->SetStatusChange(fSuspend);
	if(xBiasOn)
	{
	  G4Exception("Cannot cross-section bias a process that suspends tracks.");
	}
      }
      if(aR->GetStatusChange()!=stopAndKill )
      {
	if(xBiasOn && G4UniformRand()<XBiasSurvivalProbability() )
	{
 	  theTotalResult->SetWeightChange( XBiasSurvivalProbability() );
	  G4double newWeight = aR->GetWeightChange();
          G4DynamicParticle * aNew = new G4DynamicParticle(aT.GetDefinition(),
	                                                   aR->GetEnergyChange(),
							   aR->GetMomentumChange());
	  G4HadSecondary * theSec = new G4HadSecondary(aNew, newWeight);
	  aR->AddSecondary(theSec);
	}
	else
	{
	  G4double newWeight = aR->GetWeightChange();
	  theTotalResult->SetWeightChange(newWeight); // This is multiplicative
	  if(aR->GetEnergyChange()>-.5) theTotalResult->SetEnergyChange(aR->GetEnergyChange());
	  G4LorentzVector newDirection(aR->GetMomentumChange().unit(), 1.);
	  newDirection*=aR->GetTrafoToLab();
	  theTotalResult->SetMomentumDirectionChange(newDirection.vect());
	}
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
	if(xBiasOn) newWeight *= XBiasSecondaryWeight();
	track->SetWeight(newWeight);
	theTotalResult->AddSecondary(track);
      }
      return;
}
 /* end of file */
