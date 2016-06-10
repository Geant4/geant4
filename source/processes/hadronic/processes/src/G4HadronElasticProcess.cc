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
// $Id: G4HadronElasticProcess.cc 92396 2015-08-31 14:12:40Z gcosmo $
//
// Geant4 Hadron Elastic Scattering Process 
// 
// Created 26 July 2012 V.Ivanchenko from G4WHadronElasticProcess
//  
// Modified:
// 14-Sep-12 M.Kelsey -- Pass subType code to base ctor

#include <iostream>
#include <typeinfo>

#include "G4HadronElasticProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4Nucleus.hh"
#include "G4ProcessManager.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4HadronElasticDataSet.hh"
#include "G4ProductionCutsTable.hh"
#include "G4HadronicException.hh"
#include "G4HadronicDeprecate.hh"
#include "G4HadronicInteraction.hh"
#include "G4VCrossSectionRatio.hh"

G4HadronElasticProcess::G4HadronElasticProcess(const G4String& pName)
  : G4HadronicProcess(pName, fHadronElastic), isInitialised(false),
    fDiffraction(0), fDiffractionRatio(0)
{
  AddDataSet(new G4HadronElasticDataSet);
  lowestEnergy = 1.*keV;
}

G4HadronElasticProcess::~G4HadronElasticProcess()
{}

void G4HadronElasticProcess::ProcessDescription(std::ostream& outFile) const
{

    outFile << "G4HadronElasticProcess handles the elastic scattering of \n"
            << "hadrons by invoking the following hadronic model(s) and \n"
            << "hadronic cross section(s).\n";

}

G4VParticleChange* 
G4HadronElasticProcess::PostStepDoIt(const G4Track& track, 
				     const G4Step& /*step*/)
{
  theTotalResult->Clear();
  theTotalResult->Initialize(track);
  G4double weight = track.GetWeight();
  theTotalResult->ProposeWeight(weight);

  // For elastic scattering, _any_ result is considered an interaction
  ClearNumberOfInteractionLengthLeft();

  G4double kineticEnergy = track.GetKineticEnergy();
  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();
  const G4ParticleDefinition* part = dynParticle->GetDefinition();

  // NOTE:  Very low energy scatters were causing numerical (FPE) errors
  //        in earlier releases; these limits have not been changed since.
  if (kineticEnergy <= lowestEnergy)   return theTotalResult;

  G4Material* material = track.GetMaterial();
  G4Nucleus* targNucleus = GetTargetNucleusPointer();

  // Select element
  G4Element* elm = 0;
  try
    {
      elm = GetCrossSectionDataStore()->SampleZandA(dynParticle, material, 
						    *targNucleus);
    }
  catch(G4HadronicException & aR)
    {
      G4ExceptionDescription ed;
      aR.Report(ed);
      DumpState(track,"SampleZandA",ed); 
      ed << " PostStepDoIt failed on element selection" << G4endl;
      G4Exception("G4HadronElasticProcess::PostStepDoIt", "had003", 
		  FatalException, ed);
    }

  // Initialize the hadronic projectile from the track
  //  G4cout << "track " << track.GetDynamicParticle()->Get4Momentum()<<G4endl;
  G4HadProjectile theProj(track);
  G4HadronicInteraction* hadi = 0;
  G4HadFinalState* result = 0;

  if(fDiffraction) 
  {
    G4double ratio = 
      fDiffractionRatio->ComputeRatio(part, kineticEnergy,
				      targNucleus->GetZ_asInt(),
				      targNucleus->GetA_asInt());
    // diffraction is chosen
    if(ratio > 0.0 && G4UniformRand() < ratio) 
    {
      try
	{
//   if(part->GetParticleName() == "pi-")
// G4cout<<part->GetParticleName()<<"; "<<kineticEnergy/CLHEP::GeV<<" GeV; r = "<<ratio<<G4endl;

	  result = fDiffraction->ApplyYourself(theProj, *targNucleus);
	}
      catch(G4HadronicException aR)
	{
	  G4ExceptionDescription ed;
	  aR.Report(ed);
	  ed << "Call for " << fDiffraction->GetModelName() << G4endl;
	  ed << "Target element "<< elm->GetName()<<"  Z= " 
	     << targNucleus->GetZ_asInt() 
	     << "  A= " << targNucleus->GetA_asInt() << G4endl;
	  DumpState(track,"ApplyYourself",ed);
	  ed << " ApplyYourself failed" << G4endl;
	  G4Exception("G4HadronElasticProcess::PostStepDoIt", "had006", 
		      FatalException, ed);
	}
      // Check the result for catastrophic energy non-conservation
      result = CheckResult(theProj, *targNucleus, result);

      result->SetTrafoToLab(theProj.GetTrafoToLab());
      ClearNumberOfInteractionLengthLeft();

      FillResult(result, track);

      if (epReportLevel != 0) {
	CheckEnergyMomentumConservation(track, *targNucleus);
      }
      return theTotalResult;
    }
  }      

  // ordinary elastic scattering
  try
    {
      hadi = ChooseHadronicInteraction( theProj, *targNucleus, material, elm );
    }
  catch(G4HadronicException & aE)
    {
      G4ExceptionDescription ed;
      aE.Report(ed);
      ed << "Target element "<< elm->GetName()<<"  Z= " 
	 << targNucleus->GetZ_asInt() << "  A= " 
	 << targNucleus->GetA_asInt() << G4endl;
      DumpState(track,"ChooseHadronicInteraction",ed);
      ed << " No HadronicInteraction found out" << G4endl;
      G4Exception("G4HadronElasticProcess::PostStepDoIt", "had005", 
		  FatalException, ed);
    }

  size_t idx = track.GetMaterialCutsCouple()->GetIndex();
  G4double tcut = (*(G4ProductionCutsTable::GetProductionCutsTable()
		     ->GetEnergyCutsVector(3)))[idx];
  hadi->SetRecoilEnergyThreshold(tcut);

  if(verboseLevel>1) {
    G4cout << "G4HadronElasticProcess::PostStepDoIt for " 
	   << part->GetParticleName()
	   << " in " << material->GetName() 
	   << " Target Z= " << targNucleus->GetZ_asInt() 
	   << " A= " << targNucleus->GetA_asInt() << G4endl; 
  }

  try
    {
      result = hadi->ApplyYourself( theProj, *targNucleus);
    }
  catch(G4HadronicException aR)
    {
      G4ExceptionDescription ed;
      aR.Report(ed);
      ed << "Call for " << hadi->GetModelName() << G4endl;
      ed << "Target element "<< elm->GetName()<<"  Z= " 
	 << targNucleus->GetZ_asInt() 
	 << "  A= " << targNucleus->GetA_asInt() << G4endl;
      DumpState(track,"ApplyYourself",ed);
      ed << " ApplyYourself failed" << G4endl;
      G4Exception("G4HadronElasticProcess::PostStepDoIt", "had006", 
		  FatalException, ed);
    }

  // Check the result for catastrophic energy non-conservation
  // cannot be applied because is not guranteed that recoil 
  // nucleus is created
  // result = CheckResult(theProj, targNucleus, result);

  // directions
  G4ThreeVector indir = track.GetMomentumDirection();
  G4double phi = CLHEP::twopi*G4UniformRand();
  G4ThreeVector it(0., 0., 1.);
  G4ThreeVector outdir = result->GetMomentumChange();

  if(verboseLevel>1) {
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
  if(efinal <= lowestEnergy) {
    edep += efinal;
    efinal = 0.0;
  }

  // primary change
  theTotalResult->ProposeEnergy(efinal);

  G4TrackStatus status = track.GetTrackStatus();
  if(efinal > 0.0) {
    outdir.rotate(phi, it);
    outdir.rotateUz(indir);
    theTotalResult->ProposeMomentumDirection(outdir);
  } else {
    if(part->GetProcessManager()->GetAtRestProcessVector()->size() > 0)
         { status = fStopButAlive; }
    else { status = fStopAndKill; }
    theTotalResult->ProposeTrackStatus(status);
  }

  //G4cout << "Efinal= " << efinal << "  TrackStatus= " << status << G4endl;

  theTotalResult->SetNumberOfSecondaries(0);

  // recoil
  if(result->GetNumberOfSecondaries() > 0) {
    G4DynamicParticle* p = result->GetSecondary(0)->GetParticle();

    if(p->GetKineticEnergy() > tcut) {
      theTotalResult->SetNumberOfSecondaries(1);
      G4ThreeVector pdir = p->GetMomentumDirection();
      // G4cout << "recoil " << pdir << G4endl;
      //!! is not needed for models inheriting G4HadronElastic
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

    } else {
      edep += p->GetKineticEnergy();
      delete p;
    }
  }
  theTotalResult->ProposeLocalEnergyDeposit(edep);
  theTotalResult->ProposeNonIonizingEnergyDeposit(edep);
  result->Clear();

  return theTotalResult;
}

void 
G4HadronElasticProcess::PreparePhysicsTable(const G4ParticleDefinition& part)
{
  if(!isInitialised) {
    isInitialised = true;
    if(G4Neutron::Neutron() == &part) { lowestEnergy = 1.e-6*eV; }
  }
  G4HadronicProcess::PreparePhysicsTable(part);
}

void G4HadronElasticProcess::SetLowestEnergy(G4double val)
{
  lowestEnergy = val;
}

void 
G4HadronElasticProcess::SetLowestEnergyNeutron(G4double val)
{
  lowestEnergy = val;
  G4HadronicDeprecate("G4HadronElasticProcess::SetLowestEnergyNeutron()");
}

void G4HadronElasticProcess::SetDiffraction(G4HadronicInteraction* hi, 
					    G4VCrossSectionRatio* xsr)
{
  if(hi && xsr) {
    fDiffraction = hi;
    fDiffractionRatio = xsr;
  }
}
