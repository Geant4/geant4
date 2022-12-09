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
// Geant4 Hadron Inelastic Scattering Process 
// 
// Created  from G4MuNeutrinoNucleusProcess
//  


#include <iostream>
#include <typeinfo>

#include "G4TauNeutrinoNucleusProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4Nucleus.hh"
#include "G4ProcessManager.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4ProductionCutsTable.hh"
#include "G4HadronicException.hh"
#include "G4HadronicInteraction.hh"
#include "G4VCrossSectionRatio.hh"
#include "G4VDiscreteProcess.hh"

#include "G4TauNeutrinoNucleusTotXsc.hh"
//#include "G4NuMuNucleusCcModel.hh"
//#include "G4NuMuNucleusNcModel.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"
#include "G4DynamicParticle.hh"
#include "G4StepPoint.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4SafetyHelper.hh"
#include "G4TransportationManager.hh"

///////////////////////////////////////////////////////////////////////////////


G4TauNeutrinoNucleusProcess::G4TauNeutrinoNucleusProcess( G4String anEnvelopeName, const G4String& pName)
  : G4HadronicProcess( pName, fHadronInelastic ), isInitialised(false), fBiased(true)  // fHadronElastic???
{
  lowestEnergy = 1.*keV;
  fEnvelope  = nullptr;
  fEnvelopeName = anEnvelopeName;
  fTotXsc = nullptr; // new G4TauNeutrinoNucleusTotXsc();
  fNuNuclCcBias=1.;
  fNuNuclNcBias=1.;
  fNuNuclTotXscBias=1.;
  safetyHelper = G4TransportationManager::GetTransportationManager()->GetSafetyHelper();
  safetyHelper->InitialiseHelper();
}

G4TauNeutrinoNucleusProcess::~G4TauNeutrinoNucleusProcess()
{
  if( fTotXsc ) delete fTotXsc;
}

///////////////////////////////////////////////////////

void G4TauNeutrinoNucleusProcess::SetBiasingFactor(G4double bf)
{
  fNuNuclTotXscBias = bf;

  fTotXsc = new G4TauNeutrinoNucleusTotXsc();
  fTotXsc->SetBiasingFactor(bf);
}

///////////////////////////////////////////////////////

void G4TauNeutrinoNucleusProcess::SetBiasingFactors(G4double bfCc, G4double bfNc)
{
  fNuNuclCcBias = bfCc;
  fNuNuclNcBias = bfNc;

  fTotXsc = new G4TauNeutrinoNucleusTotXsc();
  // fTotXsc->SetBiasingFactors(bfCc, bfNc);
}

//////////////////////////////////////////////////

G4double G4TauNeutrinoNucleusProcess::
GetMeanFreePath(const G4Track &aTrack, G4double, G4ForceCondition *)
{
  //G4cout << "GetMeanFreePath " << aTrack.GetDefinition()->GetParticleName()
  //	 << " Ekin= " << aTrack.GetKineticEnergy() << G4endl;
  G4String rName = aTrack.GetStep()->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetRegion()->GetName();
  G4double totxsc(0.);

  if( rName == fEnvelopeName && fNuNuclTotXscBias > 1.)    
  {
      totxsc = fNuNuclTotXscBias*
	GetCrossSectionDataStore()->ComputeCrossSection(aTrack.GetDynamicParticle(),
						    aTrack.GetMaterial());
  }
  else
  {
      totxsc = GetCrossSectionDataStore()->ComputeCrossSection(aTrack.GetDynamicParticle(),
						    aTrack.GetMaterial());
  }
  G4double res = (totxsc>0.0) ? 1.0/totxsc : DBL_MAX;
  //G4cout << "         xsection= " << totxsc << G4endl;
  return res;
}

///////////////////////////////////////////////////

void G4TauNeutrinoNucleusProcess::ProcessDescription(std::ostream& outFile) const
{

    outFile << "G4TauNeutrinoNucleusProcess handles the inelastic scattering of \n"
            << "tau-neutrino on nucleus by invoking the following  model(s) and \n"
            << "cross section(s).\n";

}

///////////////////////////////////////////////////////////////////////

G4VParticleChange* 
G4TauNeutrinoNucleusProcess::PostStepDoIt(const G4Track& track, const G4Step& step)
{
  // track.GetVolume()->GetLogicalVolume()->GetName()
  // if( track.GetVolume()->GetLogicalVolume() != fEnvelope )
 
  G4String rName = track.GetStep()->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetRegion()->GetName();

  if( rName != fEnvelopeName ) 
  {
    if( verboseLevel > 0 )
    {
      G4cout<<"Go out from G4TauNeutrinoNucleusProcess::PostStepDoIt: wrong volume "<<G4endl;
    }
    return G4VDiscreteProcess::PostStepDoIt( track,  step );
  }
  theTotalResult->Clear();
  theTotalResult->Initialize(track);
  G4double weight = track.GetWeight();
  theTotalResult->ProposeWeight(weight);

  if( track.GetTrackStatus() != fAlive ) 
  { 
    return theTotalResult; 
  }
  // Next check for illegal track status
  //
  if (track.GetTrackStatus() != fAlive && 
      track.GetTrackStatus() != fSuspend) 
  {
    if (track.GetTrackStatus() == fStopAndKill ||
        track.GetTrackStatus() == fKillTrackAndSecondaries ||
        track.GetTrackStatus() == fPostponeToNextEvent) 
    {
      G4ExceptionDescription ed;
      ed << "G4TauNeutrinoNucleusProcess: track in unusable state - "
	 << track.GetTrackStatus() << G4endl;
      ed << "G4TauNeutrinoNucleusProcess: returning unchanged track " << G4endl;
      DumpState(track,"PostStepDoIt",ed);
      G4Exception("G4TauNeutrinoNucleusProcess::PostStepDoIt", "had004", JustWarning, ed);
    }
    // No warning for fStopButAlive which is a legal status here
    return theTotalResult;
  }
    
  // For elastic scattering, _any_ result is considered an interaction
  ClearNumberOfInteractionLengthLeft();

  G4double kineticEnergy = track.GetKineticEnergy();
  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();
  const G4ParticleDefinition* part = dynParticle->GetDefinition();
  const G4String pName = part->GetParticleName();

  // NOTE:  Very low energy scatters were causing numerical (FPE) errors
  //        in earlier releases; these limits have not been changed since.

  if ( kineticEnergy <= lowestEnergy )   return theTotalResult;

  const G4Material* material = track.GetMaterial();
  G4Nucleus* targNucleus = GetTargetNucleusPointer();

  //////////////// uniform random spread of the neutrino interaction point ////////////

  const G4StepPoint* pPostStepPoint  = step.GetPostStepPoint();
  const G4DynamicParticle* aParticle = track.GetDynamicParticle();
  G4ThreeVector      position  = pPostStepPoint->GetPosition(), newPosition=position;
  G4ParticleMomentum direction = aParticle->GetMomentumDirection();
  
  if( fNuNuclCcBias > 1.0 ||  fNuNuclNcBias > 1.0) // = true, if fBiasingfactor != 1., i.e. xsc is biased
  {
    const G4RotationMatrix* rotM = pPostStepPoint->GetTouchable()->GetRotation();
    G4ThreeVector transl = pPostStepPoint->GetTouchable()->GetTranslation();
    G4AffineTransform transform = G4AffineTransform(rotM,transl);
    transform.Invert();

    G4ThreeVector localP = transform.TransformPoint(position);
    G4ThreeVector localV = transform.TransformAxis(direction);

    G4double forward  = track.GetVolume()->GetLogicalVolume()->GetSolid()->DistanceToOut(localP,  localV);
    G4double backward = track.GetVolume()->GetLogicalVolume()->GetSolid()->DistanceToOut(localP, -localV);

    G4double distance = forward+backward;

    // G4cout<<distance/cm<<", ";

    // uniform sampling of nu-e interaction point 
    // along neutrino direction in current volume

    G4double range = -backward+G4UniformRand()*distance;

    newPosition = position + range*direction;

    safetyHelper->ReLocateWithinVolume(newPosition);

    theTotalResult->ProposePosition(newPosition); // G4Exception : GeomNav1002
  }
  G4HadProjectile theProj( track );
  G4HadronicInteraction* hadi = nullptr;
  G4HadFinalState* result = nullptr;

  G4double ccTotRatio = fTotXsc->GetCcTotRatio();

  if( G4UniformRand() < ccTotRatio )  // Cc-model
  {
    // Initialize the hadronic projectile from the track
    thePro.Initialise(track);

    if (pName == "nu_tau" ) hadi = (GetHadronicInteractionList())[0];
    else                    hadi = (GetHadronicInteractionList())[2];

    result = hadi->ApplyYourself( thePro, *targNucleus);
   
    result->SetTrafoToLab(thePro.GetTrafoToLab());

    ClearNumberOfInteractionLengthLeft();

    FillResult(result, track);
  }
  else  // Nc-model                            
  {

    if (pName == "nu_tau" )    hadi = (GetHadronicInteractionList())[1]; 
    else                      hadi = (GetHadronicInteractionList())[3];

    size_t idx = track.GetMaterialCutsCouple()->GetIndex();

    G4double tcut = (*(G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(3)))[idx];

    hadi->SetRecoilEnergyThreshold(tcut);

    if( verboseLevel > 1 ) 
    {
    G4cout << "G4TauNeutrinoNucleusProcess::PostStepDoIt for " 
	   << part->GetParticleName()
	   << " in " << material->GetName() 
	   << " Target Z= " << targNucleus->GetZ_asInt() 
	   << " A= " << targNucleus->GetA_asInt() << G4endl; 
    }
    try
    {
      result = hadi->ApplyYourself( theProj, *targNucleus);
    }
    catch(G4HadronicException & aR)
    {
      G4ExceptionDescription ed;
      aR.Report(ed);
      ed << "Call for " << hadi->GetModelName() << G4endl;
      ed << "  Z= " 
	 << targNucleus->GetZ_asInt() 
	 << "  A= " << targNucleus->GetA_asInt() << G4endl;
      DumpState(track,"ApplyYourself",ed);
      ed << " ApplyYourself failed" << G4endl;
      G4Exception("G4TauNeutrinoNucleusProcess::PostStepDoIt", "had006", 
		  FatalException, ed);
    }
    // directions

    G4ThreeVector indir = track.GetMomentumDirection();
    G4double phi = CLHEP::twopi*G4UniformRand();
    G4ThreeVector it(0., 0., 1.);
    G4ThreeVector outdir = result->GetMomentumChange();

    if(verboseLevel>1) 
    {
      G4cout << "Efin= " << result->GetEnergyChange()
	   << " de= " << result->GetLocalEnergyDeposit()
	   << " nsec= " << result->GetNumberOfSecondaries()
	   << " dir= " << outdir
	   << G4endl;
    }
    // energies 
 
    G4double edep = result->GetLocalEnergyDeposit();
    G4double efinal = result->GetEnergyChange();

    if(efinal < 0.0) { efinal = 0.0; }
    if(edep < 0.0)   { edep = 0.0; }

    // NOTE:  Very low energy scatters were causing numerical (FPE) errors
    //        in earlier releases; these limits have not been changed since.

    if(efinal <= lowestEnergy) 
    {
      edep += efinal;
      efinal = 0.0;
    }
    // primary change

    theTotalResult->ProposeEnergy(efinal);

    G4TrackStatus status = track.GetTrackStatus();

    if(efinal > 0.0) 
    {
      outdir.rotate(phi, it);
      outdir.rotateUz(indir);
      theTotalResult->ProposeMomentumDirection(outdir);
    } 
    else 
    {
      if( part->GetProcessManager()->GetAtRestProcessVector()->size() > 0)
      { 
        status = fStopButAlive; 
      }
      else 
      { 
        status = fStopAndKill; 
      }
      theTotalResult->ProposeTrackStatus(status);
    }
    //G4cout << "Efinal= " << efinal << "  TrackStatus= " << status << G4endl;

    theTotalResult->SetNumberOfSecondaries(0);

    // recoil

    if( result->GetNumberOfSecondaries() > 0 ) 
    {
      G4DynamicParticle* p = result->GetSecondary(0)->GetParticle();

      if(p->GetKineticEnergy() > tcut) 
      {
        theTotalResult->SetNumberOfSecondaries(1);
        G4ThreeVector pdir = p->GetMomentumDirection();

        // G4cout << "recoil " << pdir << G4endl;
        //!! is not needed for models inheriting G4TauNeutrinoNucleus

        pdir.rotate(phi, it);
        pdir.rotateUz(indir);

        // G4cout << "recoil rotated " << pdir << G4endl;

        p->SetMomentumDirection(pdir);

        // in elastic scattering time and weight are not changed

        G4Track* t = new G4Track(p, track.GetGlobalTime(), 
			       track.GetPosition());
        t->SetWeight(weight);
        t->SetTouchableHandle(track.GetTouchableHandle());
        theTotalResult->AddSecondary(t);
      } 
      else 
      {
        edep += p->GetKineticEnergy();
        delete p;
      }
    }
    theTotalResult->ProposeLocalEnergyDeposit(edep);
    theTotalResult->ProposeNonIonizingEnergyDeposit(edep);
    result->Clear();
  }
  return theTotalResult;
}

void 
G4TauNeutrinoNucleusProcess::PreparePhysicsTable(const G4ParticleDefinition& part)
{
  if(!isInitialised) {
    isInitialised = true;
    // if(G4Neutron::Neutron() == &part) { lowestEnergy = 1.e-6*eV; }
  }
  G4HadronicProcess::PreparePhysicsTable(part);
}

void 
G4TauNeutrinoNucleusProcess::SetLowestEnergy(G4double val)
{
  lowestEnergy = val;
}

