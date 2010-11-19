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
// 13.01.10: M.Kosov: Commented not used G4QElasticCrossSection & G4QCHIPSWorld
//

#include "G4WHadronElasticProcess.hh"
#include "globals.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4HadronElasticDataSet.hh"
#include "G4VQCrossSection.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4IsotopeVector.hh"
#include "G4Neutron.hh"
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
      return G4VDiscreteProcess::PostStepDoIt(track,step);
  } else {
    if(kineticEnergy <= lowestEnergy)
      return G4VDiscreteProcess::PostStepDoIt(track,step);
  }

  G4Material* material = track.GetMaterial();
  G4CrossSectionDataStore* store = GetCrossSectionDataStore();
  G4double xsec = store->GetCrossSection(dynParticle,material);
  if(xsec <= DBL_MIN) return G4VDiscreteProcess::PostStepDoIt(track,step);

  // Select element
  G4Element* elm = store->SampleZandA(dynParticle,material,targetNucleus);

  G4HadronicInteraction* hadi = 
    ChooseHadronicInteraction( kineticEnergy, material, elm);

  size_t idx = track.GetMaterialCutsCouple()->GetIndex();
  G4double tcut = 
    (*(G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(3)))[idx];
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
  G4ThreeVector outdir = (result->GetMomentumChange()).rotateUz(indir);
  
  if(verboseLevel>1) {
    G4cout << "Efin= " << result->GetEnergyChange()
	   << " de= " << result->GetLocalEnergyDeposit()
	   << " nsec= " << result->GetNumberOfSecondaries()
	   << " dir= " << outdir
	   << G4endl;
  }
  
  aParticleChange.ProposeEnergy(result->GetEnergyChange());
  aParticleChange.ProposeMomentumDirection(outdir);
  if(result->GetNumberOfSecondaries() > 0) {
    aParticleChange.SetNumberOfSecondaries(1);
    G4DynamicParticle* p = result->GetSecondary(0)->GetParticle();
    G4ThreeVector pdir = p->GetMomentumDirection();
    // G4cout << "recoil " << pdir << G4endl;
    pdir = pdir.rotateUz(indir);
    // G4cout << "recoil rotated " << pdir << G4endl;
    p->SetMomentumDirection(pdir);
    aParticleChange.AddSecondary(p);
  } else {
    aParticleChange.SetNumberOfSecondaries(0);
    aParticleChange.ProposeLocalEnergyDeposit(result->GetLocalEnergyDeposit());
  }
  result->Clear();

  return G4VDiscreteProcess::PostStepDoIt(track,step);
}

G4bool G4WHadronElasticProcess::
IsApplicable(const G4ParticleDefinition& aParticleType)
{
   return (&aParticleType == G4PionPlus::PionPlus() ||
           &aParticleType == G4PionMinus::PionMinus() ||
           &aParticleType == G4KaonPlus::KaonPlus() ||
           &aParticleType == G4KaonZeroShort::KaonZeroShort() ||
           &aParticleType == G4KaonZeroLong::KaonZeroLong() ||
           &aParticleType == G4KaonMinus::KaonMinus() ||
           &aParticleType == G4Proton::Proton() ||
           &aParticleType == G4AntiProton::AntiProton() ||
           &aParticleType == G4Neutron::Neutron() ||
           &aParticleType == G4AntiNeutron::AntiNeutron() ||
           &aParticleType == G4Lambda::Lambda() ||
           &aParticleType == G4AntiLambda::AntiLambda() ||
           &aParticleType == G4SigmaPlus::SigmaPlus() ||
           &aParticleType == G4SigmaZero::SigmaZero() ||
           &aParticleType == G4SigmaMinus::SigmaMinus() ||
           &aParticleType == G4AntiSigmaPlus::AntiSigmaPlus() ||
           &aParticleType == G4AntiSigmaZero::AntiSigmaZero() ||
           &aParticleType == G4AntiSigmaMinus::AntiSigmaMinus() ||
           &aParticleType == G4XiZero::XiZero() ||
           &aParticleType == G4XiMinus::XiMinus() ||
           &aParticleType == G4AntiXiZero::AntiXiZero() ||
           &aParticleType == G4AntiXiMinus::AntiXiMinus() ||
           &aParticleType == G4Deuteron::Deuteron() ||
           &aParticleType == G4Triton::Triton() ||
           &aParticleType == G4He3::He3() ||
           &aParticleType == G4Alpha::Alpha() ||
           &aParticleType == G4OmegaMinus::OmegaMinus() ||
           &aParticleType == G4AntiOmegaMinus::AntiOmegaMinus() ||
           &aParticleType == G4GenericIon::GenericIon());
}


