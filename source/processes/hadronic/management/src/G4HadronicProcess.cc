//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4HadronicProcess.cc,v 1.93 2010-12-01 02:04:39 dennis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class source file
//
// G4HadronicProcess
//
// original by H.P.Wellisch
// J.L. Chuma, TRIUMF, 10-Mar-1997
//
// Modifications:
// 05-Jul-2010 V.Ivanchenko cleanup commented lines 
// 20-Jul-2011 M.Kelsey -- null-pointer checks in DumpState()
// 24-Sep-2011 M.Kelsey -- Use envvar G4HADRONIC_RANDOM_FILE to save random
//		engine state before each model call
// 18-Oct-2011 M.Kelsey -- Handle final-state cases in conservation checks.

#include "G4Types.hh"
#include "G4HadronicProcess.hh"

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
#include "G4NucleiProperties.hh"

#include "G4HadronicException.hh"
#include "G4HadronicProcessStore.hh"

#include <typeinfo>
#include <sstream>
//#include <stdlib.h>

// File-scope variable to capture environment variable at startup

static const char* G4Hadronic_Random_File = getenv("G4HADRONIC_RANDOM_FILE");

// Initialize static variables for isotope production

G4IsoParticleChange * G4HadronicProcess::theIsoResult = 0;
G4IsoParticleChange * G4HadronicProcess::theOldIsoResult = 0;
G4bool G4HadronicProcess::isoIsEnabled = true;

void G4HadronicProcess::
EnableIsotopeProductionGlobally()  {isoIsEnabled = true;}

void G4HadronicProcess::
DisableIsotopeProductionGlobally() {isoIsEnabled = false;}

//////////////////////////////////////////////////////////////////

G4HadronicProcess::G4HadronicProcess(const G4String& processName,
                                     G4ProcessType aType)
 :G4VDiscreteProcess(processName, aType)
{
  ModelingState = 0;
  isoIsOnAnyway = -1;
  theTotalResult = new G4ParticleChange();
  theTotalResult->SetSecondaryWeightByProcess(true);
  theInteraction = 0;
  theCrossSectionDataStore = new G4CrossSectionDataStore();
  G4HadronicProcessStore::Instance()->Register(this);
  aScaleFactor = 1;
  xBiasOn = false;
  G4HadronicProcess_debug_flag = false;
  epReportLevel = 0;
  epCheckLevels.first = DBL_MAX;
  epCheckLevels.second = DBL_MAX;
  levelsSetByProcess = false;

  // Make ep checking possible via environment variables
  if ( char * ReportLevel = getenv("G4Hadronic_epReportLevel")) {
     std::stringstream sRL (ReportLevel);
     sRL >> epReportLevel;
     //-GF we now take min of process and model   levelsSetByProcess = true;
     if ( char * RelativeLevel = getenv("G4Hadronic_epCheckRelativeLevel")) {
     	std::stringstream level(RelativeLevel);
	level >> epCheckLevels.first;
     }
     if ( char * AbsoluteLevel = getenv("G4Hadronic_epCheckAbsoluteLevel")) {
     	std::stringstream level(AbsoluteLevel);
	level >> epCheckLevels.second;
     }
     //G4cout << " Checking E/p with level " << epReportLevel 
     //       << ", relative/absolute level = " << epCheckLevels.first << " / "<< epCheckLevels.second << G4endl;      
  }  
}

G4HadronicProcess::~G4HadronicProcess()
{ 
  G4HadronicProcessStore::Instance()->DeRegister(this);
  delete theTotalResult;

  std::for_each(theProductionModels.begin(),
                theProductionModels.end(), G4Delete());
 
  delete theOldIsoResult; 
  delete theIsoResult;
  delete theCrossSectionDataStore;
}

void G4HadronicProcess::RegisterMe( G4HadronicInteraction *a )
{ 
  if(!a) { return; }
  try{GetManagerPointer()->RegisterMe( a );}   
  catch(G4HadronicException & aE)
  {
    G4ExceptionDescription ed;
    ed << "Unrecoverable error in " << GetProcessName() 
       << " to register " << a->GetModelName() << G4endl;
    G4Exception("G4HadronicProcess::RegisterMe", "had001", FatalException,
		ed);
  }
  G4HadronicProcessStore::Instance()->RegisterInteraction(this, a);  
}

void G4HadronicProcess::PreparePhysicsTable(const G4ParticleDefinition& p)
{
  if(getenv("G4HadronicProcess_debug")) { 
    G4HadronicProcess_debug_flag = true;
  }
  G4HadronicProcessStore::Instance()->RegisterParticle(this, &p);
}

void G4HadronicProcess::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  theCrossSectionDataStore->BuildPhysicsTable(p);
  G4HadronicProcessStore::Instance()->PrintInfo(&p);
}

G4double G4HadronicProcess::
GetMeanFreePath(const G4Track &aTrack, G4double, G4ForceCondition *)
{ 
  try
  {
    theLastCrossSection = aScaleFactor* 
      theCrossSectionDataStore->GetCrossSection(aTrack.GetDynamicParticle(), 
						aTrack.GetMaterial());
  }
  catch(G4HadronicException aR)
  { 
    G4ExceptionDescription ed;
    DumpState(aTrack,"GetMeanFreePath",ed);
    ed << " Cross section is not available" << G4endl;
    G4Exception("G4HadronicProcess::GetMeanFreePath", "had002", FatalException,
		ed);
  } 
  G4double res = DBL_MAX;
  if( theLastCrossSection > 0.0 ) { res = 1.0/theLastCrossSection; }
  return res;
}

G4VParticleChange*
G4HadronicProcess::PostStepDoIt(const G4Track& aTrack, const G4Step&)
{
  // if primary is not Alive then do nothing
  theTotalResult->Initialize(aTrack);
  if(aTrack.GetTrackStatus() != fAlive) { return theTotalResult; }

  // Find cross section at end of step and check if <= 0
  //
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  G4Material* aMaterial = aTrack.GetMaterial();
   
  G4Element* anElement = 0;
  try
  {
     anElement = theCrossSectionDataStore->SampleZandA(aParticle, 
						       aMaterial, 
						       targetNucleus);
  }
  catch(G4HadronicException & aR)
  {
    G4ExceptionDescription ed;
    DumpState(aTrack,"SampleZandA",ed); 
    ed << " PostStepDoIt failed on element selection" << G4endl;
    G4Exception("G4HadronicProcess::PostStepDoIt", "had003", FatalException,
		ed);
  }

  if(aParticle->GetDefinition()->GetPDGCharge() != 0.0) {
    if (GetElementCrossSection(aParticle, anElement, aMaterial) <= 0.0) {
      // No interaction
      //theTotalResult->Clear();
      return theTotalResult;
    }    
  }

  // Next check for illegal track status
  //
  if (aTrack.GetTrackStatus() != fAlive && aTrack.GetTrackStatus() != fSuspend) {
    if (aTrack.GetTrackStatus() == fStopAndKill ||
        aTrack.GetTrackStatus() == fKillTrackAndSecondaries ||
        aTrack.GetTrackStatus() == fPostponeToNextEvent) {
      G4ExceptionDescription ed;
      ed << "G4HadronicProcess: track in unusable state - "
	 << aTrack.GetTrackStatus() << G4endl;
      ed << "G4HadronicProcess: returning unchanged track " << G4endl;
      DumpState(aTrack,"PostStepDoIt",ed);
      G4Exception("G4HadronicProcess::PostStepDoIt", "had004", JustWarning, ed);
    }
    // No warning for fStopButAlive which is a legal status here
    // theTotalResult->Clear();
    return theTotalResult;
  }

  // Go on to regular case
  //
  G4double originalEnergy = aParticle->GetKineticEnergy();
  G4double kineticEnergy = originalEnergy;

  // Get kinetic energy per nucleon for ions
  if(aParticle->GetParticleDefinition()->GetBaryonNumber() > 1.5) 
          kineticEnergy/=aParticle->GetParticleDefinition()->GetBaryonNumber();

  try
  {
    theInteraction = 
      ChooseHadronicInteraction( kineticEnergy, aMaterial, anElement );
  }
  catch(G4HadronicException & aE)
  {
    G4ExceptionDescription ed;
    ed << "Target element "<<anElement->GetName()<<"  Z= " 
       << targetNucleus.GetZ_asInt() << "  A= " 
       << targetNucleus.GetA_asInt() << G4endl;
    DumpState(aTrack,"ChooseHadronicInteraction",ed);
    ed << " No HadronicInteraction found out" << G4endl;
    G4Exception("G4HadronicProcess::PostStepDoIt", "had005", FatalException,
		ed);
  }

  // Initialize the hadronic projectile from the track

  G4HadProjectile thePro(aTrack);
  
  G4HadFinalState* result = 0;
  G4int reentryCount = 0;

  do
  {
    try
    {
      // Save random engine if requested for debugging
      if (G4Hadronic_Random_File) {
	CLHEP::HepRandom::saveEngineStatus(G4Hadronic_Random_File);
      }
      // Call the interaction
      result = theInteraction->ApplyYourself( thePro, targetNucleus);
      ++reentryCount;
    }
    catch(G4HadronicException aR)
    {
      G4ExceptionDescription ed;
      ed << "Call for " << theInteraction->GetModelName() << G4endl;
      ed << "Target element "<<anElement->GetName()<<"  Z= " 
	 << targetNucleus.GetZ_asInt() 
	 << "  A= " << targetNucleus.GetA_asInt() << G4endl;
      DumpState(aTrack,"ApplyYourself",ed);
      ed << " ApplyYourself failed" << G4endl;
      G4Exception("G4HadronicProcess::PostStepDoIt", "had006", FatalException,
		  ed);
    }
    if(reentryCount>100) {
      G4ExceptionDescription ed;
      ed << "Call for " << theInteraction->GetModelName() << G4endl;
      ed << "Target element "<<anElement->GetName()<<"  Z= " 
	 << targetNucleus.GetZ_asInt() 
	 << "  A= " << targetNucleus.GetA_asInt() << G4endl;
      DumpState(aTrack,"ApplyYourself",ed);
      ed << " ApplyYourself does not completed after 100 attempts" << G4endl;
      G4Exception("G4HadronicProcess::PostStepDoIt", "had006", FatalException,
		  ed);  
    }
  }
  while(!result);

  result->SetTrafoToLab(thePro.GetTrafoToLab());

  ClearNumberOfInteractionLengthLeft();
  /*
  if(isoIsOnAnyway!=-1)
  {
    if(isoIsEnabled||isoIsOnAnyway)
    {
      result = DoIsotopeCounting(result, aTrack, targetNucleus);
    }
  }
  // Put hadronic final state particles into G4ParticleChange

  FillTotalResult(result, aTrack);
  */

  // VI: new method   
  FillResult(result, aTrack);

  if (epReportLevel != 0) { 
    CheckEnergyMomentumConservation(aTrack, targetNucleus);
  }
  return theTotalResult;
}


void G4HadronicProcess::ProcessDescription(std::ostream& outFile) const
{
  outFile << "The description for this process has not been written yet.\n"; 
}


G4HadFinalState* 
G4HadronicProcess::DoIsotopeCounting(G4HadFinalState * aResult,
                                     const G4Track & aTrack,
                                     const G4Nucleus & aNucleus)
{
  // get the PC from iso-production
  delete theOldIsoResult;
  theOldIsoResult = 0;
  delete theIsoResult;
  theIsoResult = new G4IsoParticleChange;
  G4bool done = false;
  G4IsoResult * anIsoResult = 0;
  for(unsigned int i=0; i<theProductionModels.size(); i++)
  {
    anIsoResult = theProductionModels[i]->GetIsotope(aTrack, aNucleus);
    if(anIsoResult!=0)
    {
      done = true;
      break;
    }
  }

  // If no production models active, use default iso production
  if(!done) anIsoResult = ExtractResidualNucleus(aTrack, aNucleus, aResult); 

  // Add all info explicitely and add typename from model called.
  theIsoResult->SetIsotope(anIsoResult->GetIsotope());
  theIsoResult->SetProductionPosition(aTrack.GetPosition());
  theIsoResult->SetProductionTime(aTrack.GetGlobalTime());
  theIsoResult->SetParentParticle(*aTrack.GetDynamicParticle());
  theIsoResult->SetMotherNucleus(anIsoResult->GetMotherNucleus());
  theIsoResult->SetProducer(typeid(*theInteraction).name());
  
  delete anIsoResult;

  // If isotope production is enabled the GetIsotopeProductionInfo() 
  // method must be called or else a memory leak will result
  //
  // The following code will fix the memory leak, but remove the 
  // isotope information:
  //
  //  if(theIsoResult) {
  //    delete theIsoResult;
  //    theIsoResult = 0;
  //  }
  
  return aResult;
}

G4IsoResult* 
G4HadronicProcess::ExtractResidualNucleus(const G4Track&,
                                          const G4Nucleus& aNucleus,
                                          G4HadFinalState* aResult)
{
  G4double A = aNucleus.GetA_asInt();
  G4double Z = aNucleus.GetZ_asInt();
  G4double bufferA = 0;
  G4double bufferZ = 0;
  
  // loop over aResult, and decrement A, Z accordingly
  // cash the max
  for(G4int i=0; i<aResult->GetNumberOfSecondaries(); ++i)
  {
    G4HadSecondary* aSecTrack = aResult->GetSecondary(i);
    const G4ParticleDefinition* part = aSecTrack->GetParticle()->GetParticleDefinition(); 
    G4double Q = part->GetPDGCharge()/eplus;
    G4double N = part->GetBaryonNumber();
    if(bufferA < N)
    {
      bufferA = N;
      bufferZ = Q;
    }
    Z -= Q;
    A -= N;
  }
  
  // if the fragment was part of the final state, it is 
  // assumed to be the heaviest secondary.
  if(A<0.1)
  {
    A = bufferA;
    Z = bufferZ;
  }
  
  // prepare the IsoResult.

  std::ostringstream ost1;
  ost1 <<Z<<"_"<<A;
  G4String biff = ost1.str();
  G4IsoResult * theResult = new G4IsoResult(biff, aNucleus);

  return theResult;
}

G4double G4HadronicProcess::XBiasSurvivalProbability()
{
  G4double result = 0;
  G4double nLTraversed = GetTotalNumberOfInteractionLengthTraversed();
  G4double biasedProbability = 1.-std::exp(-nLTraversed);
  G4double realProbability = 1-std::exp(-nLTraversed/aScaleFactor);
  result = (biasedProbability-realProbability)/biasedProbability;
  return result;
}

G4double G4HadronicProcess::XBiasSecondaryWeight()
{
  G4double result = 0;
  G4double nLTraversed = GetTotalNumberOfInteractionLengthTraversed();
  result = 
     1./aScaleFactor*std::exp(-nLTraversed/aScaleFactor*(1-1./aScaleFactor));
  return result;
}

void 
G4HadronicProcess::FillResult(G4HadFinalState * aR, const G4Track & aT)
{
  theTotalResult->Clear();
  theTotalResult->Initialize(aT);
  theTotalResult->ProposeLocalEnergyDeposit(aR->GetLocalEnergyDeposit());  

  G4double rotation = CLHEP::twopi*G4UniformRand();
  G4ThreeVector it(0., 0., 1.);

  G4double efinal = aR->GetEnergyChange();
  if(efinal < 0.0) { efinal = 0.0; }

  // check status of primary
  if(aR->GetStatusChange() == stopAndKill) {
    theTotalResult->ProposeTrackStatus(fStopAndKill);
    theTotalResult->ProposeEnergy( 0.0 );

    // check its final energy
  } else if(0.0 == efinal) {
    theTotalResult->ProposeEnergy( 0.0 );
    if(aT.GetParticleDefinition()->GetProcessManager()
       ->GetAtRestProcessVector()->size() > 0)
         { aParticleChange.ProposeTrackStatus(fStopButAlive); }
    else { aParticleChange.ProposeTrackStatus(fStopAndKill); }

    // primary is not killed apply rotation and Lorentz transformation
  } else  {
    theTotalResult->ProposeTrackStatus(fAlive);
    G4double mass = aT.GetParticleDefinition()->GetPDGMass();
    G4double newE = efinal + mass;
    G4double newP = std::sqrt(efinal*(efinal + 2*mass));
    G4ThreeVector newPV = newP*aR->GetMomentumChange();
    G4LorentzVector newP4(newE, newPV);
    newP4.rotate(rotation, it);
    newP4 *= aR->GetTrafoToLab();
    theTotalResult->ProposeMomentumDirection(newP4.vect().unit());
    newE = newP4.e() - mass;
    if(G4HadronicProcess_debug_flag && newE <= 0.0) {
      G4ExceptionDescription ed;
      DumpState(aT,"Primary has zero energy after interaction",ed);
      G4Exception("G4HadronicProcess::FillResults", "had011", JustWarning, ed);
    }
    if(newE < 0.0) { newE = 0.0; }
    theTotalResult->ProposeEnergy( newE );
  }

  // check secondaries: apply rotation and Lorentz transformation
  G4int nSec = aR->GetNumberOfSecondaries();
  theTotalResult->SetNumberOfSecondaries(nSec);
 
  if(nSec > 0) {
    G4double time0 = aT.GetGlobalTime();
    for(G4int i=0; i<nSec; ++i)
      {
	G4LorentzVector theM = aR->GetSecondary(i)->GetParticle()->Get4Momentum();
	theM.rotate(rotation, it);
	theM *= aR->GetTrafoToLab();
	aR->GetSecondary(i)->GetParticle()->Set4Momentum(theM);
	G4double time = aR->GetSecondary(i)->GetTime();
	if(time<time0) { time = time0; }

	G4Track* track = new G4Track(aR->GetSecondary(i)->GetParticle(),
				     time,
				     aT.GetPosition());
	G4double newWeight = aT.GetWeight()*aR->GetSecondary(i)->GetWeight();
	// G4cout << "#### ParticleDebug "
	// <<GetProcessName()<<" "
	// <<aR->GetSecondary(i)->GetParticle()->GetDefinition()->GetParticleName()<<" "
	// <<aScaleFactor<<" "
	// <<XBiasSurvivalProbability()<<" "
	// <<XBiasSecondaryWeight()<<" "
	// <<aT.GetWeight()<<" "
	// <<aR->GetSecondary(i)->GetWeight()<<" "
	// <<aR->GetSecondary(i)->GetParticle()->Get4Momentum()<<" "
	// <<G4endl;
	track->SetWeight(newWeight);
	track->SetTouchableHandle(aT.GetTouchableHandle());
	theTotalResult->AddSecondary(track);
	if(G4HadronicProcess_debug_flag) {
	  G4double e = track->GetKineticEnergy();
          if(e <= 0.0) {
	    G4ExceptionDescription ed;
	    DumpState(aT,"Secondary has zero energy",ed);
            ed << "Secondary " << track->GetDefinition()->GetParticleName() 
	       << G4endl;
	    G4Exception("G4HadronicProcess::FillResults", "had011", JustWarning,ed);
	  }
	}
      }
  }

  aR->Clear();
  return;
}

void 
G4HadronicProcess::FillTotalResult(G4HadFinalState * aR, const G4Track & aT)
{
  theTotalResult->Clear();
  theTotalResult->ProposeLocalEnergyDeposit(0.);
  theTotalResult->Initialize(aT);
  theTotalResult->SetSecondaryWeightByProcess(true);
  theTotalResult->ProposeTrackStatus(fAlive);
  G4double rotation = CLHEP::twopi*G4UniformRand();
  G4ThreeVector it(0., 0., 1.);

  if(aR->GetStatusChange()==stopAndKill)
  {
    if( xBiasOn && G4UniformRand()<XBiasSurvivalProbability() )
    {
      theTotalResult->ProposeParentWeight( XBiasSurvivalProbability()*aT.GetWeight() );
    }
    else
    {
      theTotalResult->ProposeTrackStatus(fStopAndKill);
      theTotalResult->ProposeEnergy( 0.0 );
    }
  }
  else if(aR->GetStatusChange()!=stopAndKill )
  {
    if(aR->GetStatusChange()==suspend)
    {
      theTotalResult->ProposeTrackStatus(fSuspend);
      if(xBiasOn)
      {
	G4ExceptionDescription ed;
        DumpState(aT,"FillTotalResult",ed);
        G4Exception("G4HadronicProcess::FillTotalResult", "had007", FatalException,
		    ed,"Cannot cross-section bias a process that suspends tracks.");
      }
    } else if (aT.GetKineticEnergy() == 0) {
      theTotalResult->ProposeTrackStatus(fStopButAlive);
    }

    if(xBiasOn && G4UniformRand()<XBiasSurvivalProbability())
    {
      theTotalResult->ProposeParentWeight( XBiasSurvivalProbability()*aT.GetWeight() );
      G4double newWeight = aR->GetWeightChange()*aT.GetWeight();
      G4double newM=aT.GetParticleDefinition()->GetPDGMass();
      G4double newE=aR->GetEnergyChange() + newM;
      G4double newP=std::sqrt(newE*newE - newM*newM);
      G4DynamicParticle * aNew = 
      new G4DynamicParticle(aT.GetParticleDefinition(), newE, newP*aR->GetMomentumChange());
      aR->AddSecondary(G4HadSecondary(aNew, newWeight));
    }
    else
    {
      G4double newWeight = aR->GetWeightChange()*aT.GetWeight();
      theTotalResult->ProposeParentWeight(newWeight); // This is multiplicative
      if(aR->GetEnergyChange()>-.5) 
      {
        theTotalResult->ProposeEnergy(aR->GetEnergyChange());
      }
      G4LorentzVector newDirection(aR->GetMomentumChange().unit(), 1.);
      newDirection*=aR->GetTrafoToLab();
      theTotalResult->ProposeMomentumDirection(newDirection.vect());
    }
  }
  else
  {
    G4ExceptionDescription ed;
    G4cout << "Call for " << theInteraction->GetModelName() << G4endl;
    G4cout << "Target Z= " 
	   << targetNucleus.GetZ_asInt() 
	   << "  A= " << targetNucleus.GetA_asInt() << G4endl;
    DumpState(aT,"FillTotalResult",ed);
    G4Exception("G4HadronicProcess", "had008", FatalException,
    "use of unsupported track-status.");
  }

  if(GetProcessName() != "hElastic" && GetProcessName() != "HadronElastic"
     &&  theTotalResult->GetTrackStatus()==fAlive
     && aR->GetStatusChange()==isAlive)
    {
    // Use for debugging:   G4double newWeight = theTotalResult->GetParentWeight();

    G4double newKE = std::max(DBL_MIN, aR->GetEnergyChange());
    G4DynamicParticle* aNew = new G4DynamicParticle(aT.GetParticleDefinition(), 
                                                    aR->GetMomentumChange(), 
                                                    newKE);
    aR->AddSecondary(aNew);
    aR->SetStatusChange(stopAndKill);

    theTotalResult->ProposeTrackStatus(fStopAndKill);
    theTotalResult->ProposeEnergy( 0.0 );

  }
  theTotalResult->ProposeLocalEnergyDeposit(aR->GetLocalEnergyDeposit());
  theTotalResult->SetNumberOfSecondaries(aR->GetNumberOfSecondaries());

  if(aR->GetStatusChange() != stopAndKill)
  {
    G4double newM=aT.GetParticleDefinition()->GetPDGMass();
    G4double newE=aR->GetEnergyChange() + newM;
    G4double newP=std::sqrt(newE*newE - newM*newM);
    G4ThreeVector newPV = newP*aR->GetMomentumChange();
    G4LorentzVector newP4(newE, newPV);
    newP4.rotate(rotation, it);
    newP4*=aR->GetTrafoToLab();
    theTotalResult->ProposeMomentumDirection(newP4.vect().unit());
  }

  for(G4int i=0; i<aR->GetNumberOfSecondaries(); ++i)
  {
    G4LorentzVector theM = aR->GetSecondary(i)->GetParticle()->Get4Momentum();
    theM.rotate(rotation, it);
    theM*=aR->GetTrafoToLab();
    aR->GetSecondary(i)->GetParticle()->Set4Momentum(theM);
    G4double time = aR->GetSecondary(i)->GetTime();
    if(time<0) time = aT.GetGlobalTime();

    G4Track* track = new G4Track(aR->GetSecondary(i)->GetParticle(),
				 time,
				 aT.GetPosition());

    G4double newWeight = aT.GetWeight()*aR->GetSecondary(i)->GetWeight();
    if(xBiasOn) { newWeight *= XBiasSecondaryWeight(); }
    // G4cout << "#### ParticleDebug "
    // <<GetProcessName()<<" "
    // <<aR->GetSecondary(i)->GetParticle()->GetDefinition()->GetParticleName()<<" "
    // <<aScaleFactor<<" "
    // <<XBiasSurvivalProbability()<<" "
    // <<XBiasSecondaryWeight()<<" "
    // <<aT.GetWeight()<<" "
    // <<aR->GetSecondary(i)->GetWeight()<<" "
    // <<aR->GetSecondary(i)->GetParticle()->Get4Momentum()<<" "
    // <<G4endl;
    track->SetWeight(newWeight);
    track->SetTouchableHandle(aT.GetTouchableHandle());
    theTotalResult->AddSecondary(track);
  }

  aR->Clear();
  return;
}

G4IsoParticleChange* G4HadronicProcess::GetIsotopeProductionInfo() 
{ 
  G4IsoParticleChange * anIsoResult = theIsoResult;
  if(theIsoResult) theOldIsoResult = theIsoResult;
  theIsoResult = 0;
  return anIsoResult;
}


void G4HadronicProcess::BiasCrossSectionByFactor(G4double aScale) 
{
  xBiasOn = true;
  aScaleFactor = aScale;
  G4String it = GetProcessName(); 
  if( (it != "PhotonInelastic") && 
      (it != "ElectroNuclear") && 
      (it != "PositronNuclear") )
    {
      G4ExceptionDescription ed;
      G4Exception("G4HadronicProcess::BiasCrossSectionByFactor", "had009", FatalException, ed,
		  "Cross-section biasing available only for gamma and electro nuclear reactions.");
    }
  if(aScale<100)
    {
      G4ExceptionDescription ed;
      G4Exception("G4HadronicProcess::BiasCrossSectionByFactor", "had010", JustWarning,ed,
		  "Cross-section bias readjusted to be above safe limit. New value is 100");
      aScaleFactor = 100.;
    }
}

void 
G4HadronicProcess::CheckEnergyMomentumConservation(const G4Track& aTrack,
                                                   const G4Nucleus& aNucleus)
{
  G4double targetMass = 
    G4NucleiProperties::GetNuclearMass(aNucleus.GetA_asInt(),aNucleus.GetZ_asInt());
  G4LorentzVector projectile4mom = aTrack.GetDynamicParticle()->Get4Momentum();
  G4LorentzVector target4mom(0, 0, 0, targetMass);
  G4LorentzVector initial4mom = projectile4mom + target4mom;

  // Compute final-state momentum for scattering and "do nothing" results
  G4LorentzVector final4mom;
  if (theTotalResult->GetTrackStatus() == fStopAndKill) {
    G4Track* sec;
    G4int nSec = theTotalResult->GetNumberOfSecondaries();
    for (G4int i = 0; i < nSec; i++) {
      sec = theTotalResult->GetSecondary(i);
      final4mom += sec->GetDynamicParticle()->Get4Momentum();
    }
  } else {	// Interaction didn't complete, returned "do nothing" state
    G4Track temp(aTrack);
    temp.SetMomentumDirection(*theTotalResult->GetMomentumDirection());
    temp.SetKineticEnergy(theTotalResult->GetEnergy());
    final4mom = temp.GetDynamicParticle()->Get4Momentum() + target4mom;
  }

  // Get level-checking information (used to cut-off relative checks)
  G4String processName = GetProcessName();
  G4HadronicInteraction* theModel = GetHadronicInteraction();
  G4String modelName("none");
  if (theModel) modelName = theModel->GetModelName();
  std::pair<G4double, G4double> checkLevels = epCheckLevels;
  if (!levelsSetByProcess) {
    if (theModel) checkLevels = theModel->GetEnergyMomentumCheckLevels();
    checkLevels.first= std::min(checkLevels.first,  epCheckLevels.first);
    checkLevels.second=std::min(checkLevels.second, epCheckLevels.second);
  }

  // Compute absolute total-energy difference, and relative kinetic-energy
  G4bool checkRelative = (aTrack.GetKineticEnergy() > checkLevels.second);

  G4LorentzVector diff = initial4mom - final4mom;
  G4double absolute = diff.e();
  G4double relative = checkRelative ? absolute/aTrack.GetKineticEnergy() : 0.;

  // Evaluate relative and absolute conservation
  G4bool relPass = false;
  G4String relResult = "fail";
  if (std::abs(relative) < checkLevels.first) {
    relPass = true;
    relResult = checkRelative ? "pass" : "N/A";
  }

  G4bool absPass = false;
  G4String absResult = "fail";
  if (std::abs(absolute) < checkLevels.second) {
    absPass = true;
    absResult = "pass";
  }

  std::stringstream Myout;
  // Options for level of reporting detail:
  //  0. off
  //  1. report only when E/p not conserved
  //  2. report regardless of E/p conservation
  //  3. report only when E/p not conserved, with model names, process names, and limits 
  //  4. report regardless of E/p conservation, with model names, process names, and limits
  //  negative -1.., as above, but send output to stderr

  if(std::abs(epReportLevel) == 4) {
    Myout << " Process: " << processName << " , Model: " <<  modelName << G4endl; 
    Myout << " relative limit " << checkLevels.first << " relative value = "
           << relative << " " << relResult << G4endl;
    Myout << " absolute limit (MeV) " << checkLevels.second/MeV << " absolute value (MeV) = "
           << absolute/MeV << " " << absResult << G4endl;

  } else if(std::abs(epReportLevel) == 3) {
    if (!absPass || !relPass) {
      Myout << " Process: " << processName << " , Model: " <<  modelName << G4endl;
      Myout << " Primary: " << aTrack.GetParticleDefinition()->GetParticleName()
            << " (" << aTrack.GetParticleDefinition()->GetPDGEncoding() << "),"
            << " E= " <<  aTrack.GetDynamicParticle()->Get4Momentum().e()
	    << ", target nucleus (" << aNucleus.GetZ_asInt() << "," 
	    << aNucleus.GetA_asInt() << ")" << G4endl;
      Myout << " relative limit " << checkLevels.first << " relative value = "
             << relative << " " << relResult << G4endl;  
      Myout << " absolute limit (MeV) " << checkLevels.second/MeV << " absolute value (MeV) = "
             << absolute/MeV << " " << absResult << G4endl;
    }

  } else if(std::abs(epReportLevel) == 2) {
    Myout << " relative value = " << relative << " " << relPass
             << " absolute value (MeV) = " << absolute/MeV << " " << absPass << G4endl;

  } else if(std::abs(epReportLevel) == 1) {
    if (!absPass || !relPass) {
      Myout << " relative value = " << relative << " " << relPass
             << " absolute value (MeV) = " << absolute/MeV << " " << absPass << G4endl;
    }
  }
  
  if (epReportLevel > 0)      G4cout << Myout.str();
  else if (epReportLevel < 0) G4cerr << Myout.str();
}


void G4HadronicProcess::DumpState(const G4Track& aTrack, 
				  const G4String& method,
				  G4ExceptionDescription& ed)
{
  ed << "Unrecoverable error in the method " << method << " of " 
     << GetProcessName() << G4endl;
  ed << "TrackID= "<< aTrack.GetTrackID() << "  ParentID= " 
     << aTrack.GetParentID()
     << "  " << aTrack.GetParticleDefinition()->GetParticleName() 
     << G4endl;
  ed << "Ekin(GeV)= " << aTrack.GetKineticEnergy()/CLHEP::GeV 
     << ";  direction= " << aTrack.GetMomentumDirection() << G4endl; 
  ed << "Position(mm)= " << aTrack.GetPosition()/CLHEP::mm << ";";

  if (aTrack.GetMaterial()) { 
    ed << "  material " << aTrack.GetMaterial()->GetName();
  }
  ed << G4endl;

  if (aTrack.GetVolume()) {
    ed << "PhysicalVolume  <" << aTrack.GetVolume()->GetName() 
       << ">" << G4endl;
  }
} 

/* end of file */
