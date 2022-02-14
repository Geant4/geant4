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
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronHElasticPhysics
//
// Author: 23 November 2006 V. Ivanchenko
//
// Modified:
// 21.03.07 (V.Ivanchenko) Use G4BGGNucleonElasticXS and G4BGGPionElasticXS; 
//                         Reduce thresholds for HE and Q-models to zero
// 03.06.2010 V.Ivanchenko cleanup constructors and ConstructProcess method
// 06.06.2014 A.Ribon      Use the current best elastic models.
//
//----------------------------------------------------------------------------
//

#include "G4HadronHElasticPhysics.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4Neutron.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronElastic.hh"
#include "G4ChipsElasticModel.hh"
#include "G4AntiNuclElastic.hh"
#include "G4DiffuseElastic.hh"
#include "G4NuclNuclDiffuseElastic.hh"

#include "G4CrossSectionElastic.hh"
#include "G4BGGNucleonElasticXS.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4NeutronElasticXS.hh"
#include "G4ChipsKaonMinusElasticXS.hh"
#include "G4ChipsKaonPlusElasticXS.hh"
#include "G4ChipsKaonZeroElasticXS.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4LMsdGenerator.hh"
#include "G4DiffElasticRatio.hh"

#include "G4HadronicParameters.hh"
#include "G4HadronicBuilder.hh"
#include "G4HadParticles.hh"
#include "G4HadProcesses.hh"
#include "G4PhysListUtil.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY( G4HadronHElasticPhysics );

G4HadronHElasticPhysics::G4HadronHElasticPhysics(G4int ver, G4bool diffraction)
  : G4HadronElasticPhysics(ver, "hElastic_BEST"), 
    fDiffraction(diffraction) 
{
  if (ver > 1) { 
    G4cout << "### G4HadronHElasticPhysics: " << GetPhysicsName() 
	   << "  low-mass diffraction: " << fDiffraction << G4endl; 
  }
}

G4HadronHElasticPhysics::~G4HadronHElasticPhysics() {}

void G4HadronHElasticPhysics::ConstructProcess() {
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();
  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  const G4double elimitDiffuse = 10.0;
  const G4double elimitAntiNuc = 100.0*MeV;
  const G4double delta = 0.1*MeV;
  G4double emax = std::max(param->GetMaxEnergy(), elimitAntiNuc+delta);
  if (param->GetVerboseLevel() > 1 ) {
    G4cout << "### HadronHElasticPhysics::ConstructProcess: lower energy limit for DiffuseElastic : " 
	   << elimitDiffuse/GeV << " GeV" << G4endl
           << "                                             transition energy for anti-nuclei : " 
	   << elimitAntiNuc/GeV << " GeV" << G4endl;
  }
  G4HadronElastic* lhep0 = new G4HadronElastic();
  G4HadronElastic* lhep1 = new G4HadronElastic();
  G4HadronElastic* lhep2 = new G4HadronElastic();
  lhep0->SetMaxEnergy(emax);
  lhep1->SetMaxEnergy(elimitDiffuse+delta);
  lhep2->SetMaxEnergy(elimitAntiNuc+delta);

  G4AntiNuclElastic* anuc = new G4AntiNuclElastic();
  anuc->SetMinEnergy( elimitAntiNuc );
  anuc->SetMaxEnergy(emax);

  auto anucxs = G4HadProcesses::ElasticXS("AntiAGlauber");
  auto xsNN = G4HadProcesses::ElasticXS("Glauber-Gribov Nucl-nucl");

  G4LMsdGenerator* diffGen = nullptr;
  G4DiffElasticRatio* diffRatio = nullptr;
  if( fDiffraction ) {
    diffGen = new G4LMsdGenerator("LMsdDiffraction");
    diffRatio = new G4DiffElasticRatio();
  }

  // Use Chips elastic model only for the hydrogen element and above an energy threshold
  G4ChipsElasticModel* chipsH = new G4ChipsElasticModel();
  chipsH->SetMinEnergy( elimitDiffuse );
  const G4ElementTable* theElementTable = G4Element::GetElementTable();
  for ( size_t i_ele = 0; i_ele < theElementTable->size(); i_ele++ ) {
    G4Element* element = (*theElementTable)[ i_ele ];
    if ( element->GetZ() > 1.0 ) chipsH->DeActivateFor( element );
  }

  // p
  G4ParticleDefinition* particle = G4Proton::Proton();
  G4HadronElasticProcess* hel = new G4HadronElasticProcess();
  hel->AddDataSet(new G4BGGNucleonElasticXS(particle));
  G4DiffuseElastic* protonDiffuseElastic = new G4DiffuseElastic();
  protonDiffuseElastic->SetMinEnergy( elimitDiffuse );
  hel->RegisterMe( chipsH );  // Use Chips only for Hydrogen element
  hel->RegisterMe( protonDiffuseElastic );
  hel->RegisterMe( lhep1 );
  if( fDiffraction) hel->SetDiffraction(diffGen, diffRatio);
  if( useFactorXS ) hel->MultiplyCrossSectionBy( param->XSFactorNucleonElastic() );
  ph->RegisterProcess(hel, particle);

  // n
  particle = G4Neutron::Neutron();
  hel = new G4HadronElasticProcess();
  hel->AddDataSet(new G4NeutronElasticXS());
  G4DiffuseElastic* neutronDiffuseElastic = new G4DiffuseElastic();
  neutronDiffuseElastic->SetMinEnergy( elimitDiffuse );
  hel->RegisterMe( chipsH );  // Use Chips only for Hydrogen element
  hel->RegisterMe( neutronDiffuseElastic );
  hel->RegisterMe( lhep1 );
  if( fDiffraction) hel->SetDiffraction(diffGen, diffRatio);
  if( useFactorXS ) hel->MultiplyCrossSectionBy( param->XSFactorNucleonElastic() );
  ph->RegisterProcess(hel, particle);

  // pi+
  particle = G4PionPlus::PionPlus();
  hel = new G4HadronElasticProcess();
  hel->AddDataSet(new G4BGGPionElasticXS(particle));
  G4DiffuseElastic* dElastic = new G4DiffuseElastic();
  dElastic->SetMinEnergy( elimitDiffuse );
  hel->RegisterMe( chipsH );  // Use Chips only for Hydrogen element
  hel->RegisterMe( dElastic );
  hel->RegisterMe( lhep1 );
  if( fDiffraction) hel->SetDiffraction(diffGen, diffRatio);
  if( useFactorXS ) hel->MultiplyCrossSectionBy( param->XSFactorPionElastic() );
  ph->RegisterProcess(hel, particle);

  // pi-
  particle = G4PionMinus::PionMinus();
  hel = new G4HadronElasticProcess();
  hel->AddDataSet(new G4BGGPionElasticXS(particle));
  dElastic = new G4DiffuseElastic();
  dElastic->SetMinEnergy( elimitDiffuse );
  hel->RegisterMe( chipsH );  // Use Chips only for Hydrogen element
  hel->RegisterMe( dElastic );
  hel->RegisterMe( lhep1 );
  if( fDiffraction) hel->SetDiffraction(diffGen, diffRatio);
  if( useFactorXS ) hel->MultiplyCrossSectionBy( param->XSFactorPionElastic() );
  ph->RegisterProcess(hel, particle);

  // kaons
  G4HadronicBuilder::BuildElastic( G4HadParticles::GetKaons() );

  // d, t, He3, alpha
  for( auto & pdg : G4HadParticles::GetLightIons() ) {
    particle = table->FindParticle( pdg );
    if ( particle == nullptr ) { continue; }

    hel = new G4HadronElasticProcess();
    hel->AddDataSet(xsNN);
    hel->RegisterMe(lhep0);
    if( useFactorXS ) hel->MultiplyCrossSectionBy( param->XSFactorHadronElastic() );
    ph->RegisterProcess(hel, particle);
  }

  // high energy particles
  if( emax > param->EnergyThresholdForHeavyHadrons() ) {

    // pbar, nbar, anti light ions
    for( auto & pdg : G4HadParticles::GetLightAntiIons() ) {
      particle = table->FindParticle( pdg );
      if ( particle == nullptr ) { continue; }

      hel = new G4HadronElasticProcess();
      hel->RegisterMe(lhep2);
      hel->RegisterMe(anuc);
      hel->AddDataSet(anucxs);
      if( useFactorXS ) hel->MultiplyCrossSectionBy( param->XSFactorHadronElastic() );
      ph->RegisterProcess(hel, particle);
    }

    // hyperons
    G4HadronicBuilder::BuildElastic( G4HadParticles::GetHyperons() );
    G4HadronicBuilder::BuildElastic( G4HadParticles::GetAntiHyperons() );

    // b-, c- baryons and mesons
    if( G4HadronicParameters::Instance()->EnableBCParticles() ) {
      G4HadronicBuilder::BuildElastic( G4HadParticles::GetBCHadrons() );
    }
  }
}
