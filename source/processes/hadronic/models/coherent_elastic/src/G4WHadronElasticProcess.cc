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
// $Id: G4WHadronElasticProcess.cc,v 1.5 2010-11-19 18:50:03 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Geant4 Hadron Elastic Scattering Process 
// 
// Created 21 April 2006 V.Ivanchenko
//  
// Modified:
// 24.04.06 V.Ivanchenko add neutron scattering on hydrogen from CHIPS
// 07.06.06 V.Ivanchenko fix problem of rotation of final state
// 25.07.06 V.Ivanchenko add 19 MeV low energy for CHIPS
// 26.09.06 V.Ivanchenko add lowestEnergy
// 20.10.06 V.Ivanchenko initialise lowestEnergy=0 for neitrals, eV for charged
// 23.01.07 V.Ivanchenko add cross section interfaces with Z and A
// 02.05.07 V.Ivanchenko add He3
// 13.01.10 M.Kosov: Commented not used G4QElasticCrossSection & G4QCHIPSWorld
// 24.02.11 V.Ivanchenko use particle name in IfApplicable, 
//                       added anti particles for light ions
// 
//

#include "G4WHadronElasticProcess.hh"
#include "globals.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4HadronElasticDataSet.hh"
#include "G4ProductionCutsTable.hh"
 
G4WHadronElasticProcess::G4WHadronElasticProcess(const G4String& pName)
  : G4HadronicProcess(pName) 
{
  SetProcessSubType(fHadronElastic);
  AddDataSet(new G4HadronElasticDataSet);
  theNeutron  = G4Neutron::Neutron();
  lowestEnergy = 1.*keV;
  lowestEnergyNeutron = 1.e-6*eV;
}

G4WHadronElasticProcess::~G4WHadronElasticProcess()
{
}

G4VParticleChange* G4WHadronElasticProcess::PostStepDoIt(
				  const G4Track& track, 
				  const G4Step& step)
{
  aParticleChange.Initialize(track);
  G4double kineticEnergy = track.GetKineticEnergy();
  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();
  const G4ParticleDefinition* part = dynParticle->GetDefinition();

  // protection against numerical problems
  if(part == theNeutron) {
    if(kineticEnergy <= lowestEnergyNeutron) 
      { return G4VDiscreteProcess::PostStepDoIt(track,step); }
  } else if(kineticEnergy <= lowestEnergy)
      { return G4VDiscreteProcess::PostStepDoIt(track,step); }

  G4Material* material = track.GetMaterial();
  G4CrossSectionDataStore* store = GetCrossSectionDataStore();
  G4double xsec = store->GetCrossSection(dynParticle,material);
  if(xsec <= 0.0) { return G4VDiscreteProcess::PostStepDoIt(track,step); }

  // Select element
  G4Element* elm = store->SampleZandA(dynParticle,material,targetNucleus);

  G4HadronicInteraction* hadi = 
    ChooseHadronicInteraction( kineticEnergy, material, elm);

  size_t idx = track.GetMaterialCutsCouple()->GetIndex();
  G4double tcut = (*(G4ProductionCutsTable::GetProductionCutsTable()
		     ->GetEnergyCutsVector(3)))[idx];
  hadi->SetRecoilEnergyThreshold(tcut);

  // Initialize the hadronic projectile from the track
  //  G4cout << "track " << track.GetDynamicParticle()->Get4Momentum()<<G4endl;
  G4HadProjectile thePro(track);
  if(verboseLevel>1) {
    G4cout << "G4WHadronElasticProcess::PostStepDoIt for " 
	   << part->GetParticleName()
	   << " in " << material->GetName() 
	   << " Target Z= " << targetNucleus.GetZ_asInt() 
	   << " A= " << targetNucleus.GetA_asInt() << G4endl; 
  }
  G4HadFinalState* result = hadi->ApplyYourself(thePro, targetNucleus);
  G4ThreeVector indir = track.GetMomentumDirection();

  //!! is not needed for models inheriting G4HadronElastic
  G4ThreeVector outdir = (result->GetMomentumChange()).rotateUz(indir);
  
  if(verboseLevel>1) {
    G4cout << "Efin= " << result->GetEnergyChange()
	   << " de= " << result->GetLocalEnergyDeposit()
	   << " nsec= " << result->GetNumberOfSecondaries()
	   << " dir= " << outdir
	   << G4endl;
  }

  // primary change  
  aParticleChange.ProposeEnergy(result->GetEnergyChange());
  aParticleChange.ProposeMomentumDirection(outdir);

  // recoil
  if(result->GetNumberOfSecondaries() > 0) {
    aParticleChange.SetNumberOfSecondaries(1);
    G4DynamicParticle* p = result->GetSecondary(0)->GetParticle();
    G4ThreeVector pdir = p->GetMomentumDirection();
    // G4cout << "recoil " << pdir << G4endl;
    //!! is not needed for models inheriting G4HadronElastic
    pdir = pdir.rotateUz(indir);
    // G4cout << "recoil rotated " << pdir << G4endl;
    p->SetMomentumDirection(pdir);
    aParticleChange.AddSecondary(p);
  } else {
    aParticleChange.SetNumberOfSecondaries(0);
    aParticleChange.ProposeLocalEnergyDeposit(result->GetLocalEnergyDeposit());
    aParticleChange.ProposeNonIonizingEnergyDeposit(result->GetLocalEnergyDeposit());
  }
  result->Clear();

  return G4VDiscreteProcess::PostStepDoIt(track,step);
}

G4bool G4WHadronElasticProcess::
IsApplicable(const G4ParticleDefinition& aParticleType)
{
  G4String particleName = aParticleType.GetParticleName();
  G4bool res = false;

  if (particleName == "pi+" ||
      particleName == "pi-" ||
      particleName == "kaon+" ||
      particleName == "kaon-" ||
      particleName == "neutron" ||
      particleName == "proton" ) { res = true; }
  else if (
	   particleName == "anti_He3" ||
	   particleName == "anti_alpha" ||
	   particleName == "anti_deuteron" ||
	   particleName == "anti_lambda" ||
	   particleName == "anti_omega-" ||
	   particleName == "anti_proton" ||
	   particleName == "anti_sigma0" ||
	   particleName == "anti_sigma+" ||
	   particleName == "anti_sigma-" ||
	   particleName == "anti_triton" ||
	   particleName == "anti_xi0" ||
	   particleName == "anti_xi-" ||
	   particleName == "deuteron" ||
	   particleName == "kaon0L" ||
	   particleName == "kaon0S" ||
	   particleName == "lambda" ||
	   particleName == "omega-" ||
	   particleName == "sigma0" ||
	   particleName == "sigma+" ||
	   particleName == "sigma-" ||
	   particleName == "tau+" ||
	   particleName == "tau-" ||
	   particleName == "triton" ||
	   particleName == "xi0" ||
	   particleName == "xi-" ) { res = true; }
  return res;  
}


