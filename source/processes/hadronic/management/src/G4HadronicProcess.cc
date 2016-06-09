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
//

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
    G4cout << "Unrecoverable error in " << GetProcessName() 
	   << " to register " << a->GetModelName() << G4endl;
    G4Exception("G4HadronicProcess", "007", FatalException,
    "Could not register G4HadronicInteraction");
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
    DumpState(aTrack,"GetMeanFreePath");
    G4Exception("G4HadronicProcess", "007", FatalException,
    "G4HadronicProcess::GetMeanFreePath failed");
  } 
  G4double res = DBL_MAX;
  if( theLastCrossSection > 0.0 ) { res = 1.0/theLastCrossSection; }
  return res;
}

G4double G4HadronicProcess::
GetMicroscopicCrossSection(const G4DynamicParticle *aParticle, 
			   const G4Element *anElement, 
			   G4double aTemp )
{
  G4double x =
    theCrossSectionDataStore->GetCrossSection(aParticle, anElement, aTemp);
  if(x < 0.0) { x = 0.0; }
  return x;
}


G4VParticleChange*
G4HadronicProcess::PostStepDoIt(const G4Track& aTrack, const G4Step&)
{
  // Find cross section at end of step and check if <= 0
  //
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  G4Material* aMaterial = aTrack.GetMaterial();
  G4double aTemp = aMaterial->GetTemperature();
   
  G4Element* anElement = 0;
  try
  {
     anElement = theCrossSectionDataStore->SampleZandA(aParticle, 
						      aMaterial, 
						      targetNucleus);
  }
  catch(G4HadronicException & aR)
  {
    DumpState(aTrack,"SampleZandA");
    G4Exception("G4HadronicProcess", "007", FatalException,
    "PostStepDoIt failed on element selection.");
  }

  // Next check for illegal track status
  //
  if (aTrack.GetTrackStatus() != fAlive && aTrack.GetTrackStatus() != fSuspend) {
    if (aTrack.GetTrackStatus() == fStopAndKill ||
        aTrack.GetTrackStatus() == fKillTrackAndSecondaries ||
        aTrack.GetTrackStatus() == fPostponeToNextEvent) {
      G4cout << "G4HadronicProcess: track in unusable state - "
             << aTrack.GetTrackStatus() << G4endl;
      G4cout << "G4HadronicProcess: returning unchanged track " << G4endl;
      DumpState(aTrack,"PostStepDoIt");
      G4Exception("G4HadronicProcess", "001", JustWarning, "bailing out");
    }
    // No warning for fStopButAlive which is a legal status here
    theTotalResult->Clear();
    theTotalResult->Initialize(aTrack);
    return theTotalResult;
  }

  if (GetMicroscopicCrossSection(aParticle, anElement, aTemp) <= 0.0) {
    // No interaction
    theTotalResult->Clear();
    theTotalResult->Initialize(aTrack);
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
    G4cout << "Target element "<<anElement->GetName()<<"  Z= " 
	   << targetNucleus.GetZ() << "  A= " << targetNucleus.GetN() << G4endl;
    DumpState(aTrack,"ChooseHadronicInteraction");
    G4Exception("G4HadronicProcess", "007", FatalException,
    "ChooseHadronicInteraction failed.");
  }

  // Initialize the hadronic projectile from the track

  G4HadProjectile thePro(aTrack);
  
  G4HadFinalState* result = 0;
  G4int reentryCount = 0;

  do
  {
    try
    {
      // Call the interaction
      result = theInteraction->ApplyYourself( thePro, targetNucleus);
      ++reentryCount;
    }
    catch(G4HadronicException aR)
    {
      G4cout << "Call for " << theInteraction->GetModelName() << G4endl;
      G4cout << "Target element "<<anElement->GetName()<<"  Z= " 
	     << targetNucleus.GetZ() << "  A= " << targetNucleus.GetN() << G4endl;
      DumpState(aTrack,"ApplyYourself");
      G4Exception("G4HadronicProcess", "007", FatalException,
      "PostStepDoIt failed.");
    }
    if(reentryCount>100) {
      G4cout << "Call for " << theInteraction->GetModelName() << G4endl;
      G4cout << "Target element "<<anElement->GetName()<<"  Z= " 
	     << targetNucleus.GetZ() << "  A= " << targetNucleus.GetN() << G4endl;
      DumpState(aTrack,"ApplyYourself");
      G4Exception("G4HadronicProcess", "007", FatalException,
		  "Reentering ApplyYourself too often - PostStepDoIt failed.");  
    }
  }
  while(!result);

  result->SetTrafoToLab(thePro.GetTrafoToLab());

  ClearNumberOfInteractionLengthLeft();
  if(isoIsOnAnyway!=-1)
  {
    if(isoIsEnabled||isoIsOnAnyway)
    {
      result = DoIsotopeCounting(result, aTrack, targetNucleus);
    }
  }
  
  // Put hadronic final state particles into G4ParticleChange

  FillTotalResult(result, aTrack);

  if (epReportLevel > 0) { 
    CheckEnergyMomentumConservation(aTrack, targetNucleus);
  }
  return theTotalResult;
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
  G4double A = aNucleus.GetN();
  G4double Z = aNucleus.GetZ();
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
G4HadronicProcess::FillTotalResult(G4HadFinalState * aR, const G4Track & aT)
{
  //  G4Nancheck go_wild;
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
        DumpState(aT,"FillTotalResult");
        G4Exception("G4HadronicProcess", "007", FatalException,
        "Cannot cross-section bias a process that suspends tracks.");
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
      G4HadSecondary * theSec = new G4HadSecondary(aNew, newWeight);
      aR->AddSecondary(theSec);
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
    G4cout << "Call for " << theInteraction->GetModelName() << G4endl;
    G4cout << "Target Z= " 
	   << targetNucleus.GetZ() << "  A= " << targetNucleus.GetN() << G4endl;
    DumpState(aT,"FillTotalResult");
    G4Exception("G4HadronicProcess", "007", FatalException,
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
    G4HadSecondary* theSec = new G4HadSecondary(aNew, 1.0);
    aR->AddSecondary(theSec);
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
    //static G4double pinelcount=0;
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
      G4Exception("G4HadronicProcess::BiasCrossSectionByFactor", "007", FatalException,
		  "Cross-section biasing available only for gamma and electro nuclear reactions.");
    }
  if(aScale<100)
    {
      G4Exception("G4HadronicProcess::BiasCrossSectionByFactor", "001", JustWarning,
		  "Cross-section bias readjusted to be above safe limit. New value is 100");
      aScaleFactor = 100.;
    }
}

void 
G4HadronicProcess::CheckEnergyMomentumConservation(const G4Track& aTrack,
                                                   const G4Nucleus& aNucleus)
{
  G4double targetMass = G4NucleiProperties::GetNuclearMass(aNucleus.GetN(), aNucleus.GetZ());
  G4LorentzVector projectile4mom = aTrack.GetDynamicParticle()->Get4Momentum();
  G4LorentzVector target4mom(0, 0, 0, targetMass);
  G4LorentzVector initial4mom = projectile4mom + target4mom;

  G4Track* sec;
  G4LorentzVector final4mom;
  G4int nSec = theTotalResult->GetNumberOfSecondaries();
  for (G4int i = 0; i < nSec; i++) {
    sec = theTotalResult->GetSecondary(i);
    final4mom += sec->GetDynamicParticle()->Get4Momentum();
  }

  G4LorentzVector diff = initial4mom - final4mom;
  G4double absolute = diff.e();
  G4double relative = absolute/aTrack.GetKineticEnergy();

  G4String processName = GetProcessName();
  G4HadronicInteraction* theModel = GetHadronicInteraction();
  G4String modelName("none");
  if (theModel) modelName = theModel->GetModelName();

  std::pair<G4double, G4double> checkLevels = epCheckLevels;;
  if (!levelsSetByProcess) {
    if (theModel) checkLevels = theModel->GetEnergyMomentumCheckLevels();
  }

  G4bool relPass = false;
  G4String relResult = "fail";
  if (std::abs(relative) < checkLevels.first) {
    relPass = true;
    relResult = "pass";
  }

  G4bool absPass = false;
  G4String absResult = "fail";
  if (std::abs(absolute) < checkLevels.second) {
    absPass = true;
    absResult = "pass";
  }

  // Options for level of reporting detail:
  //  0. off
  //  1. report only when E/p not conserved
  //  2. report regardless of E/p conservation
  //  3. report only when E/p not conserved, with model names, process names, and limits 
  //  4. report regardless of E/p conservation, with model names, process names, and limits

  if(epReportLevel == 4) {
    G4cout << " Process: " << processName << " , Model: " <<  modelName << G4endl; 
    G4cout << " relative limit " << checkLevels.first << " relative value = "
           << relative << " " << relResult << G4endl;
    G4cout << " absolute limit " << checkLevels.second << " absolute value = "
           << absolute << " " << absResult << G4endl;

  } else if(epReportLevel == 3) {
    if (!absPass || !relPass) {
      G4cout << " Process: " << processName << " , Model: " <<  modelName << G4endl; 
      G4cout << " relative limit " << checkLevels.first << " relative value = "
             << relative << " " << relResult << G4endl;  
      G4cout << " absolute limit " << checkLevels.second << " absolute value = "
             << absolute << " " << absResult << G4endl;
    }

  } else if(epReportLevel == 2) {
    G4cout << " relative value = " << relative << " " << relPass
           << " absolute value = " << absolute << " " << absPass << G4endl;

  } else if(epReportLevel == 1) {
    if (!absPass || !relPass) {
      G4cout << " relative value = " << relative << " " << relPass
             << " absolute value = " << absolute << " " << absPass << G4endl;
    }
  }
}


void G4HadronicProcess::DumpState(const G4Track& aTrack, const G4String& method)
{
  G4cout << "Unrecoverable error in the method " << method << " of " 
	 << GetProcessName() << G4endl;
  G4cout << "TrackID= "<< aTrack.GetTrackID() << "  ParentID= " << aTrack.GetParentID()
	 << "  " << aTrack.GetParticleDefinition()->GetParticleName() << G4endl;
  G4cout << "Ekin(GeV)= " << aTrack.GetKineticEnergy()/CLHEP::GeV 
	 << ";  direction= " << aTrack.GetMomentumDirection() << G4endl; 
  G4cout << "Position(mm)= " << aTrack.GetPosition()/CLHEP::mm 
	 << ";  material " << aTrack.GetMaterial()->GetName() << G4endl; 
  G4cout << "PhysicalVolume  <" << aTrack.GetVolume()->GetName() << ">" << G4endl;
} 

/* end of file */
