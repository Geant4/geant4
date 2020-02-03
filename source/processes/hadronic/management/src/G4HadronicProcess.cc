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
// 14-Mar-2012 G.Folger -- enhance checks for conservation of energy, etc.
// 28-Jul-2012 M.Maire  -- add function GetTargetDefinition() 
// 14-Sep-2012 Inherit from RestDiscrete, use subtype code (now in ctor) to
//		configure base-class
// 28-Sep-2012 Restore inheritance from G4VDiscreteProcess, remove enable-flag
//		changing, remove warning message from original ctor.
// 21-Aug-2019 V.Ivanchenko leave try/catch only for ApplyYourself(..), cleanup 

#include "G4HadronicProcess.hh"

#include "G4Types.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadProjectile.hh"
#include "G4ElementVector.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Element.hh"
#include "G4ParticleChange.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessManager.hh"
#include "G4NucleiProperties.hh"

#include "G4HadronicException.hh"
#include "G4HadronicProcessStore.hh"
#include "G4VCrossSectionDataSet.hh"

#include "G4NistManager.hh"
#include "G4PhysicsModelCatalog.hh"
#include "G4VLeadingParticleBiasing.hh"
#include "G4Exp.hh"

#include <typeinfo>
#include <sstream>
#include <iostream>

// File-scope variable to capture environment variable at startup

static const char* G4Hadronic_Random_File = std::getenv("G4HADRONIC_RANDOM_FILE");

//////////////////////////////////////////////////////////////////

G4HadronicProcess::G4HadronicProcess(const G4String& processName,
                                     G4ProcessType procType)
 : G4VDiscreteProcess(processName, procType)
{
  SetProcessSubType(fHadronInelastic);	// Default unless subclass changes
  InitialiseLocal();
}

//////////////////////////////////////////////////////////////////

G4HadronicProcess::G4HadronicProcess(const G4String& processName,
                                     G4HadronicProcessType aHadSubType)
 : G4VDiscreteProcess(processName, fHadronic)
{
  SetProcessSubType(aHadSubType);
  InitialiseLocal();
}

G4HadronicProcess::~G4HadronicProcess()
{
  theProcessStore->DeRegister(this);
  delete theTotalResult;
  delete theCrossSectionDataStore;
}

void G4HadronicProcess::InitialiseLocal() {  
  theTotalResult = new G4ParticleChange();
  theTotalResult->SetSecondaryWeightByProcess(true);
  theInteraction = nullptr;
  theCrossSectionDataStore = new G4CrossSectionDataStore();
  theProcessStore = G4HadronicProcessStore::Instance();
  theProcessStore->Register(this);
  theInitialNumberOfInteractionLength = 0.0;
  aScaleFactor = 1.0;
  fWeight = 1.0;
  xBiasOn = false;
  nMatWarn = 0;
  useIntegralXS = true;
  theLastCrossSection = 0.0;
  nICelectrons = 0;
  idxIC = -1;
  G4HadronicProcess_debug_flag = false;
  GetEnergyMomentumCheckEnvvars();
}

void G4HadronicProcess::GetEnergyMomentumCheckEnvvars() {
  levelsSetByProcess = false;

  epReportLevel = std::getenv("G4Hadronic_epReportLevel") ?
    std::strtol(std::getenv("G4Hadronic_epReportLevel"),0,10) : 0;

  epCheckLevels.first = std::getenv("G4Hadronic_epCheckRelativeLevel") ?
    std::strtod(std::getenv("G4Hadronic_epCheckRelativeLevel"),0) : DBL_MAX;

  epCheckLevels.second = std::getenv("G4Hadronic_epCheckAbsoluteLevel") ?
    std::strtod(std::getenv("G4Hadronic_epCheckAbsoluteLevel"),0) : DBL_MAX;
}

void G4HadronicProcess::RegisterMe( G4HadronicInteraction *a )
{
  if(!a) { return; }
  try{ theEnergyRangeManager.RegisterMe( a ); }
  catch(G4HadronicException & aE)
  {
    G4ExceptionDescription ed;
    aE.Report(ed);
    ed << "Unrecoverable error in " << GetProcessName()
       << " to register " << a->GetModelName() << G4endl;
    G4Exception("G4HadronicProcess::RegisterMe", "had001", FatalException,
		ed);
  }
  G4HadronicProcessStore::Instance()->RegisterInteraction(this, a);
}

G4double 
G4HadronicProcess::GetElementCrossSection(const G4DynamicParticle * part,
					  const G4Element * elm, 
					  const G4Material* mat)
{
  if(!mat) 
  {
    ++nMatWarn;
    static const G4int nmax = 5;
    if(nMatWarn < nmax) {
      G4ExceptionDescription ed;
      ed << "Cannot compute Element x-section for " << GetProcessName()
	 << " because no material defined \n"
	 << " Please, specify material pointer or define simple material"
	 << " for Z= " << elm->GetZasInt();
      G4Exception("G4HadronicProcess::GetElementCrossSection", "had066", 
		  JustWarning, ed);
    }
  }
  return
    std::max(theCrossSectionDataStore->GetCrossSection(part, elm, mat),0.0);
}

void G4HadronicProcess::PreparePhysicsTable(const G4ParticleDefinition& p)
{
  if(std::getenv("G4HadronicProcess_debug")) {
    G4HadronicProcess_debug_flag = true;
  }
  theProcessStore->RegisterParticle(this, &p);
}

void G4HadronicProcess::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  theCrossSectionDataStore->BuildPhysicsTable(p);
  theEnergyRangeManager.BuildPhysicsTable(p);
  G4HadronicProcessStore::Instance()->PrintInfo(&p);
}

G4double G4HadronicProcess::
GetMeanFreePath(const G4Track &aTrack, G4double, G4ForceCondition *)
{
  //G4cout << "GetMeanFreePath " << aTrack.GetDefinition()->GetParticleName()
  //	 << " Ekin= " << aTrack.GetKineticEnergy() << G4endl;
  theLastCrossSection = aScaleFactor*theCrossSectionDataStore
     ->ComputeCrossSection(aTrack.GetDynamicParticle(),aTrack.GetMaterial());
  G4double res = (theLastCrossSection>0.0) ? 1.0/theLastCrossSection : DBL_MAX;
  //G4cout << "         xsection= " << theLastCrossSection << G4endl;
  return res;
}

G4VParticleChange*
G4HadronicProcess::PostStepDoIt(const G4Track& aTrack, const G4Step&)
{
  //G4cout << "PostStepDoIt " << aTrack.GetDefinition()->GetParticleName()
  //	 << " Ekin= " << aTrack.GetKineticEnergy() << G4endl;
  // if primary is not Alive then do nothing
  theTotalResult->Clear();
  theTotalResult->Initialize(aTrack);
  fWeight = aTrack.GetWeight();
  theTotalResult->ProposeWeight(fWeight);
  if(aTrack.GetTrackStatus() != fAlive) { return theTotalResult; }

  // Find cross section at end of step and check if <= 0
  //
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4Material* aMaterial = aTrack.GetMaterial();

  // check only for charged particles
  if(aParticle->GetDefinition()->GetPDGCharge() != 0.0) {
    G4double xs = aScaleFactor*
      theCrossSectionDataStore->ComputeCrossSection(aParticle,aMaterial);
    if(xs <= 0.0 || xs < theLastCrossSection*G4UniformRand()) {
      // No interaction
      return theTotalResult;
    }    
  }

  const G4Element* anElement = 
    theCrossSectionDataStore->SampleZandA(aParticle,aMaterial,targetNucleus);

  // Next check for illegal track status
  //
  if (aTrack.GetTrackStatus() != fAlive && 
      aTrack.GetTrackStatus() != fSuspend) {
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
    return theTotalResult;
  }

  // Initialize the hadronic projectile from the track
  thePro.Initialise(aTrack);

  theInteraction = ChooseHadronicInteraction(thePro, targetNucleus, 
                                             aMaterial, anElement);
  if(!theInteraction) {
    G4ExceptionDescription ed;
    ed << "Target element "<<anElement->GetName()<<"  Z= "
       << targetNucleus.GetZ_asInt() << "  A= "
       << targetNucleus.GetA_asInt() << G4endl;
    DumpState(aTrack,"ChooseHadronicInteraction",ed);
    ed << " No HadronicInteraction found out" << G4endl;
    G4Exception("G4HadronicProcess::PostStepDoIt", "had005", FatalException, ed);
    return theTotalResult;
  }

  G4HadFinalState* result = nullptr;
  G4int reentryCount = 0;
  /*
  G4cout << "### " << aParticle->GetDefinition()->GetParticleName() 
	 << "  Ekin(MeV)= " << aParticle->GetKineticEnergy()
	 << "  Z= " << targetNucleus.GetZ_asInt() 
	 << "  A= " << targetNucleus.GetA_asInt()
	 << "  by " << theInteraction->GetModelName() 
	 << G4endl;
  */
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
    catch(G4HadronicException & aR)
    {
      G4ExceptionDescription ed;
      aR.Report(ed);
      ed << "Call for " << theInteraction->GetModelName() << G4endl;
      ed << "Target element "<<anElement->GetName()<<"  Z= "
	 << targetNucleus.GetZ_asInt()
	 << "  A= " << targetNucleus.GetA_asInt() << G4endl;
      DumpState(aTrack,"ApplyYourself",ed);
      ed << " ApplyYourself failed" << G4endl;
      G4Exception("G4HadronicProcess::PostStepDoIt", "had006", FatalException,
		  ed);
    }

    // Check the result for catastrophic energy non-conservation
    result = CheckResult(thePro, targetNucleus, result);

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
  while(!result);  /* Loop checking, 30-Oct-2015, G.Folger */

  // Check whether kaon0 or anti_kaon0 are present between the secondaries: 
  // if this is the case, transform them into either kaon0S or kaon0L,
  // with equal, 50% probability, keeping their dynamical masses (and
  // the other kinematical properties). 
  // When this happens - very rarely - a "JustWarning" exception is thrown.
  G4int nSec = result->GetNumberOfSecondaries();
  if ( nSec > 0 ) {
    for ( G4int i = 0; i < nSec; ++i ) {
      G4DynamicParticle* dynamicParticle = result->GetSecondary(i)->GetParticle();
      const G4ParticleDefinition* particleDefinition = 
        dynamicParticle->GetParticleDefinition();
      if ( particleDefinition == G4KaonZero::Definition() || 
           particleDefinition == G4AntiKaonZero::Definition() ) {
        G4ParticleDefinition* newPart;
	if( G4UniformRand() > 0.5 ) { newPart = G4KaonZeroShort::Definition(); }
	else { newPart = G4KaonZeroLong::Definition(); }
        dynamicParticle->SetDefinition( newPart );
        G4ExceptionDescription ed;
        ed << " Hadronic model " << theInteraction->GetModelName() << G4endl;
        ed << " created " << particleDefinition->GetParticleName() << G4endl;
        ed << " -> forced to be " << newPart->GetParticleName() << G4endl;
        G4Exception( "G4HadronicProcess::PostStepDoIt", "had007", JustWarning, ed );
      }
    }
  }

  result->SetTrafoToLab(thePro.GetTrafoToLab());

  ClearNumberOfInteractionLengthLeft();

  FillResult(result, aTrack);

  if (epReportLevel != 0) {
    CheckEnergyMomentumConservation(aTrack, targetNucleus);
  }
  //G4cout << "PostStepDoIt done " << G4endl;
  return theTotalResult;
}

void G4HadronicProcess::ProcessDescription(std::ostream& outFile) const
{
  outFile << "The description for this process has not been written yet.\n";
}


G4double G4HadronicProcess::XBiasSurvivalProbability()
{
  G4double nLTraversed = GetTotalNumberOfInteractionLengthTraversed();
  G4double biasedProbability = 1.-G4Exp(-nLTraversed);
  G4double realProbability = 1-G4Exp(-nLTraversed/aScaleFactor);
  G4double result = (biasedProbability-realProbability)/biasedProbability;
  return result;
}

G4double G4HadronicProcess::XBiasSecondaryWeight()
{
  G4double nLTraversed = GetTotalNumberOfInteractionLengthTraversed();
  G4double result =
     1./aScaleFactor*G4Exp(-nLTraversed/aScaleFactor*(1-1./aScaleFactor));
  return result;
}

void
G4HadronicProcess::FillResult(G4HadFinalState * aR, const G4Track & aT)
{
  theTotalResult->ProposeLocalEnergyDeposit(aR->GetLocalEnergyDeposit());
  const G4ThreeVector& dir = aT.GetMomentumDirection();

  G4double efinal = std::max(aR->GetEnergyChange(), 0.0);

  // check status of primary
  if(aR->GetStatusChange() == stopAndKill) {
    theTotalResult->ProposeTrackStatus(fStopAndKill);
    theTotalResult->ProposeEnergy( 0.0 );

    // check its final energy
  } else if(0.0 == efinal) {
    theTotalResult->ProposeEnergy( 0.0 );
    if(aT.GetParticleDefinition()->GetProcessManager()
       ->GetAtRestProcessVector()->size() > 0)
         { theTotalResult->ProposeTrackStatus(fStopButAlive); }
    else { theTotalResult->ProposeTrackStatus(fStopAndKill); }

    // primary is not killed apply rotation and Lorentz transformation
  } else  {
    theTotalResult->ProposeTrackStatus(fAlive);
    G4ThreeVector newDir = aR->GetMomentumChange();
    newDir.rotateUz(dir);
    theTotalResult->ProposeMomentumDirection(newDir);
    theTotalResult->ProposeEnergy(efinal);
  }
  //G4cout << "FillResult: Efinal= " << efinal << " status= " 
  //	 << theTotalResult->GetTrackStatus() 
  //	 << "  fKill= " << fStopAndKill << G4endl;

  // check secondaries 
  nICelectrons = 0;
  if(idxIC == -1) { 
    G4int idx = G4PhysicsModelCatalog::GetIndex("e-InternalConvertion");
    idxIC =  -1 == idx ? -2 : idx; 
  } 
  G4int nSec = aR->GetNumberOfSecondaries();
  theTotalResult->SetNumberOfSecondaries(nSec);
  G4double time0 = aT.GetGlobalTime();

  for (G4int i = 0; i < nSec; ++i) {
    G4DynamicParticle* dynParticle = aR->GetSecondary(i)->GetParticle();

    // apply rotation
    G4ThreeVector newDir = dynParticle->GetMomentumDirection();
    newDir.rotateUz(dir);
    dynParticle->SetMomentumDirection(newDir);

    // check if secondary is on the mass shell
    const G4ParticleDefinition* part = dynParticle->GetDefinition();
    G4double mass = part->GetPDGMass();
    G4double dmass= dynParticle->GetMass();
    const G4double delta_mass_lim = 1.0*CLHEP::keV;
    const G4double delta_ekin = 0.001*CLHEP::eV;
    if(std::abs(dmass - mass) > delta_mass_lim) {
      G4double e = std::max(dynParticle->GetKineticEnergy() + dmass - mass, delta_ekin);
      if(G4HadronicProcess_debug_flag) {
	G4ExceptionDescription ed;
	ed << "TrackID= "<< aT.GetTrackID()
	   << "  " << aT.GetParticleDefinition()->GetParticleName()
	   << " Target Z= " << targetNucleus.GetZ_asInt() << "  A= "
	   << targetNucleus.GetA_asInt() 
	   << " Ekin(GeV)= " << aT.GetKineticEnergy()/CLHEP::GeV
	   << "\n Secondary is out of mass shell: " << part->GetParticleName()
	   << "  EkinNew(MeV)= " << e  
	   << " DeltaMass(MeV)= " << dmass - mass << G4endl;
	G4Exception("G4HadronicProcess::FillResults", "had012", JustWarning, ed);
      }
      dynParticle->SetKineticEnergy(e);
      dynParticle->SetMass(mass);               
    }
    G4int idxModel = aR->GetSecondary(i)->GetCreatorModelType(); 
    if(idxIC == idxModel) { ++nICelectrons; }
      
    // time of interaction starts from zero + global time
    G4double time = std::max(aR->GetSecondary(i)->GetTime(), 0.0) + time0;

    G4Track* track = new G4Track(dynParticle, time, aT.GetPosition());
    track->SetCreatorModelIndex(idxModel);
    G4double newWeight = fWeight*aR->GetSecondary(i)->GetWeight();
    track->SetWeight(newWeight);
    track->SetTouchableHandle(aT.GetTouchableHandle());
    theTotalResult->AddSecondary(track);
    if (G4HadronicProcess_debug_flag) {
      G4double e = dynParticle->GetKineticEnergy();
      if (e == 0.0) {
	G4ExceptionDescription ed;
	DumpState(aT,"Secondary has zero energy",ed);
	ed << "Secondary " << part->GetParticleName()
	   << G4endl;
	G4Exception("G4HadronicProcess::FillResults", "had011", 
		      JustWarning,ed);
      }
    }
  }
  aR->Clear();
}

void G4HadronicProcess::MultiplyCrossSectionBy(G4double factor)
{
  BiasCrossSectionByFactor(factor);
}

void G4HadronicProcess::BiasCrossSectionByFactor(G4double aScale)
{
  if (aScale <= 0.0) {
    G4ExceptionDescription ed;
    ed << " Wrong biasing factor " << aScale << " for " << GetProcessName();
    G4Exception("G4HadronicProcess::BiasCrossSectionByFactor", "had010", 
                JustWarning, ed, "Cross-section bias is ignored");
  } else {
    xBiasOn = true;
    aScaleFactor = aScale;
  }
}

G4HadFinalState* G4HadronicProcess::CheckResult(const G4HadProjectile & aPro,
						const G4Nucleus &aNucleus, 
						G4HadFinalState * result)
{
  // check for catastrophic energy non-conservation
  // to re-sample the interaction

  G4HadronicInteraction * theModel = GetHadronicInteraction();
  G4double nuclearMass(0);
  if (theModel) {

    // Compute final-state total energy
    G4double finalE(0.);
    G4int nSec = result->GetNumberOfSecondaries();

    nuclearMass = G4NucleiProperties::GetNuclearMass(aNucleus.GetA_asInt(),
						     aNucleus.GetZ_asInt());
    if (result->GetStatusChange() != stopAndKill) {
      // Interaction didn't complete, returned "do nothing" state 
      // and reset nucleus or the primary survived the interaction 
      // (e.g. electro-nuclear ) => keep  nucleus
      finalE=result->GetLocalEnergyDeposit() +
             aPro.GetDefinition()->GetPDGMass() + result->GetEnergyChange();
      if( nSec == 0 ){
         // Since there are no secondaries, there is no recoil nucleus.
         // To check energy balance we must neglect the initial nucleus too.
        nuclearMass=0.0;
      }
    }
    for (G4int i = 0; i < nSec; i++) {
      G4DynamicParticle *pdyn=result->GetSecondary(i)->GetParticle();
      finalE += pdyn->GetTotalEnergy();
      G4double mass_pdg=pdyn->GetDefinition()->GetPDGMass();
      G4double mass_dyn=pdyn->GetMass();
      if ( std::abs(mass_pdg - mass_dyn) > 0.1*mass_pdg + 1.*MeV ) {
        // If it is shortlived, then a difference less than 3 times the width is acceptable
        if ( pdyn->GetDefinition()->IsShortLived()  &&
             std::abs(mass_pdg - mass_dyn) < 3.0*pdyn->GetDefinition()->GetPDGWidth() ) {
          continue;
        }
	result->Clear();
	result = nullptr;
	G4ExceptionDescription desc;
	desc << "Warning: Secondary with off-shell dynamic mass detected:  " 
	     << G4endl
	     << " " << pdyn->GetDefinition()->GetParticleName()
	     << ", PDG mass: " << mass_pdg << ", dynamic mass: "
	     << mass_dyn << G4endl
	     << (epReportLevel<0 ? "abort the event" 
		 : "re-sample the interaction") << G4endl
	     << " Process / Model: " <<  GetProcessName()<< " / "
	     << theModel->GetModelName() << G4endl
	     << " Primary: " << aPro.GetDefinition()->GetParticleName()
	     << " (" << aPro.GetDefinition()->GetPDGEncoding() << "), "
	     << " E= " <<  aPro.Get4Momentum().e()
	     << ", target nucleus (" << aNucleus.GetZ_asInt() << ", "
	     << aNucleus.GetA_asInt() << ")" << G4endl;
	G4Exception("G4HadronicProcess:CheckResult()", "had012",
		    epReportLevel<0 ? EventMustBeAborted : JustWarning,desc);
	// must return here.....
	return result;
      }
    }
    G4double deltaE= nuclearMass +  aPro.GetTotalEnergy() -  finalE;

    std::pair<G4double, G4double> checkLevels = 
      theModel->GetFatalEnergyCheckLevels();	// (relative, absolute)
    if (std::abs(deltaE) > checkLevels.second && 
        std::abs(deltaE) > checkLevels.first*aPro.GetKineticEnergy()){
      // do not delete result, this is a pointer to a data member;
      result->Clear();
      result = nullptr;
      G4ExceptionDescription desc;
      desc << "Warning: Bad energy non-conservation detected, will "
	   << (epReportLevel<0 ? "abort the event" 
	       :  "re-sample the interaction") << G4endl
	   << " Process / Model: " <<  GetProcessName()<< " / " 
	   << theModel->GetModelName() << G4endl
	   << " Primary: " << aPro.GetDefinition()->GetParticleName()
	   << " (" << aPro.GetDefinition()->GetPDGEncoding() << "), "
	   << " E= " <<  aPro.Get4Momentum().e()
	   << ", target nucleus (" << aNucleus.GetZ_asInt() << ", "
	   << aNucleus.GetA_asInt() << ")" << G4endl
	   << " E(initial - final) = " << deltaE << " MeV." << G4endl;
      G4Exception("G4HadronicProcess:CheckResult()", "had012", 
		  epReportLevel<0 ? EventMustBeAborted : JustWarning,desc);
    }
  }
  return result;
}

void
G4HadronicProcess::CheckEnergyMomentumConservation(const G4Track& aTrack,
                                                   const G4Nucleus& aNucleus)
{
  G4int target_A=aNucleus.GetA_asInt();
  G4int target_Z=aNucleus.GetZ_asInt();
  G4double targetMass = G4NucleiProperties::GetNuclearMass(target_A,target_Z);
  G4LorentzVector target4mom(0, 0, 0, targetMass 
			     + nICelectrons*CLHEP::electron_mass_c2);

  G4LorentzVector projectile4mom = aTrack.GetDynamicParticle()->Get4Momentum();
  G4int track_A = aTrack.GetDefinition()->GetBaryonNumber();
  G4int track_Z = G4lrint(aTrack.GetDefinition()->GetPDGCharge());

  G4int initial_A = target_A + track_A;
  G4int initial_Z = target_Z + track_Z - nICelectrons;

  G4LorentzVector initial4mom = projectile4mom + target4mom;

  // Compute final-state momentum for scattering and "do nothing" results
  G4LorentzVector final4mom;
  G4int final_A(0), final_Z(0);

  G4int nSec = theTotalResult->GetNumberOfSecondaries();
  if (theTotalResult->GetTrackStatus() != fStopAndKill) {  // If it is Alive
    // Either interaction didn't complete, returned "do nothing" state
    //  or    the primary survived the interaction (e.g. electro-nucleus )

    // Interaction didn't complete, returned "do nothing" state
    //   - or suppressed recoil  (e.g. Neutron elastic )
    final4mom = initial4mom;
    final_A = initial_A;
    final_Z = initial_Z;
    if (nSec > 0 && aTrack.GetDynamicParticle()) {
      // The primary remains in final state (e.g. electro-nucleus )
      G4Track temp(aTrack);

      // Use the final energy / momentum
      temp.SetMomentumDirection(*theTotalResult->GetMomentumDirection());
      temp.SetKineticEnergy(theTotalResult->GetEnergy());
      final4mom = temp.GetDynamicParticle()->Get4Momentum();
      final_A = track_A;
      final_Z = track_Z;
      // Expect that the target nucleus will have interacted,
      //  and its products, including recoil, will be included in secondaries.
    }
  }
  if( nSec > 0 ) {
    G4Track* sec;

    for (G4int i = 0; i < nSec; i++) {
      sec = theTotalResult->GetSecondary(i);
      final4mom += sec->GetDynamicParticle()->Get4Momentum();
      final_A += sec->GetDefinition()->GetBaryonNumber();
      final_Z += G4lrint(sec->GetDefinition()->GetPDGCharge());
    }
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

  G4double absolute_mom = diff.vect().mag();
  G4double relative_mom = checkRelative ? absolute_mom/aTrack.GetMomentum().mag() : 0.;

  // Evaluate relative and absolute conservation
  G4bool relPass = true;
  G4String relResult = "pass";
  if (  std::abs(relative) > checkLevels.first
	 || std::abs(relative_mom) > checkLevels.first) {
    relPass = false;
    relResult = checkRelative ? "fail" : "N/A";
  }

  G4bool absPass = true;
  G4String absResult = "pass";
  if (   std::abs(absolute) > checkLevels.second
      || std::abs(absolute_mom) > checkLevels.second ) {
    absPass = false ;
    absResult = "fail";
  }

  G4bool chargePass = true;
  G4String chargeResult = "pass";
  if (   (initial_A-final_A)!=0
      || (initial_Z-final_Z)!=0 ) {
    chargePass = checkLevels.second < DBL_MAX ? false : true;
    chargeResult = "fail";
   }

  G4bool conservationPass = (relPass || absPass) && chargePass;

  std::stringstream Myout;
  G4bool Myout_notempty(false);
  // Options for level of reporting detail:
  //  0. off
  //  1. report only when E/p not conserved
  //  2. report regardless of E/p conservation
  //  3. report only when E/p not conserved, with model names, process names, and limits
  //  4. report regardless of E/p conservation, with model names, process names, and limits
  //  negative -1.., as above, but send output to stderr

  if(   std::abs(epReportLevel) == 4
	||	( std::abs(epReportLevel) == 3 && ! conservationPass ) ){
      Myout << " Process: " << processName << " , Model: " <<  modelName << G4endl;
      Myout << " Primary: " << aTrack.GetParticleDefinition()->GetParticleName()
            << " (" << aTrack.GetParticleDefinition()->GetPDGEncoding() << "),"
            << " E= " <<  aTrack.GetDynamicParticle()->Get4Momentum().e()
	    << ", target nucleus (" << aNucleus.GetZ_asInt() << ","
	    << aNucleus.GetA_asInt() << ")" << G4endl;
      Myout_notempty=true;
  }
  if (  std::abs(epReportLevel) == 4
	 || std::abs(epReportLevel) == 2
	 || ! conservationPass ){

      Myout << "   "<< relResult  <<" relative, limit " << checkLevels.first << ", values E/T(0) = "
             << relative << " p/p(0)= " << relative_mom  << G4endl;
      Myout << "   "<< absResult << " absolute, limit (MeV) " << checkLevels.second/MeV << ", values E / p (MeV) = "
             << absolute/MeV << " / " << absolute_mom/MeV << " 3mom: " << (diff.vect())*1./MeV <<  G4endl;
      Myout << "   "<< chargeResult << " charge/baryon number balance " << (initial_Z-final_Z) << " / " << (initial_A-final_A) << " "<<  G4endl;
      Myout_notempty=true;

  }
  Myout.flush();
  if ( Myout_notempty ) {
     if (epReportLevel > 0)      G4cout << Myout.str()<< G4endl;
     else if (epReportLevel < 0) G4cerr << Myout.str()<< G4endl;
  }
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

void G4HadronicProcess::DumpPhysicsTable(const G4ParticleDefinition& p)
{ 
  theCrossSectionDataStore->DumpPhysicsTable(p); 
}

void G4HadronicProcess::AddDataSet(G4VCrossSectionDataSet * aDataSet)
{ 
  theCrossSectionDataStore->AddDataSet(aDataSet);
}

std::vector<G4HadronicInteraction*>& 
G4HadronicProcess::GetHadronicInteractionList()
{ 
  return theEnergyRangeManager.GetHadronicInteractionList(); 
}

G4HadronicInteraction*
G4HadronicProcess::GetHadronicModel(const G4String& modelName)
{ 
  std::vector<G4HadronicInteraction*>& list
        = theEnergyRangeManager.GetHadronicInteractionList();
  for (size_t li=0; li<list.size(); li++) {
    if (list[li]->GetModelName() == modelName) return list[li];
  }
  return nullptr;
}
