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
// $Id: G4HadronHElasticPhysics.cc 106721 2017-10-20 09:46:54Z gcosmo $
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
// CHIPS for sampling scattering for p and n
// Glauber model for samplimg of high energy pi+- (E > 1GeV)
// LHEP sampling model for the other particle
// BBG cross sections for p, n and pi+- 
// LHEP cross sections for other particles

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
#include "G4ChipsHyperonElasticXS.hh"
#include "G4ChipsAntiBaryonElasticXS.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4LMsdGenerator.hh"
#include "G4DiffElasticRatio.hh"
#include "G4AutoDelete.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY( G4HadronHElasticPhysics );

G4ThreadLocal G4DiffElasticRatio* G4HadronHElasticPhysics::diffRatio = 0;

G4HadronHElasticPhysics::G4HadronHElasticPhysics( G4int ver, G4bool diffraction)
  : G4VPhysicsConstructor( "hElastic_BEST" ), verbose( ver ), 
    fDiffraction(diffraction) 
{
  if ( verbose > 1 ) { 
    G4cout << "### G4HadronHElasticPhysics: " << GetPhysicsName() 
	   << "  low-mass diffraction: " << fDiffraction << G4endl; 
  }
}

G4HadronHElasticPhysics::~G4HadronHElasticPhysics() {}


void G4HadronHElasticPhysics::ConstructParticle() {
  // G4cout << "G4HadronElasticPhysics::ConstructParticle" << G4endl;
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void G4HadronHElasticPhysics::ConstructProcess() {

  const G4double elimitDiffuse = 0.0;
  const G4double elimitAntiNuc = 100.0*MeV;
  const G4double delta = 0.1*MeV;

  if ( verbose > 1 ) {
    G4cout << "### HadronHElasticPhysics::ConstructProcess: lower energy limit for DiffuseElastic : " 
	   << elimitDiffuse/GeV << " GeV" << G4endl
           << "                                             transition energy for anti-nuclei : " 
	   << elimitAntiNuc/GeV << " GeV" << G4endl;
  }

  G4AntiNuclElastic* anuc = new G4AntiNuclElastic();
  anuc->SetMinEnergy( elimitAntiNuc );
  G4CrossSectionElastic* anucxs = 
    new G4CrossSectionElastic( anuc->GetComponentCrossSection() );

  G4HadronElastic* lhep = new G4HadronElastic();
  lhep->SetMaxEnergy( elimitAntiNuc + delta );

  // Three instances of Chips elastic model: one used everywhere,
  // one used below a energy threshold, and one used only for the
  // hydrogen element.
  G4ChipsElasticModel* chips1 = new G4ChipsElasticModel();
  G4ChipsElasticModel* chips2 = new G4ChipsElasticModel();
  chips2->SetMaxEnergy( elimitAntiNuc + delta );
  G4ChipsElasticModel* chipsH = new G4ChipsElasticModel();
  const G4ElementTable* theElementTable = G4Element::GetElementTable();
  for ( size_t i_ele = 0; i_ele < theElementTable->size(); i_ele++ ) {
    G4Element* element = (*theElementTable)[ i_ele ];
    if ( element->GetZ() > 1.0 ) chipsH->DeActivateFor( element );
  }

  G4NuclNuclDiffuseElastic* diffuseNuclNuclElastic = new G4NuclNuclDiffuseElastic();
  diffuseNuclNuclElastic->SetMinEnergy( elimitDiffuse );

  G4VCrossSectionDataSet* theComponentGGHadronNucleusData = 
    new G4CrossSectionElastic( new G4ComponentGGHadronNucleusXsc );

  G4VCrossSectionDataSet* theComponentGGNuclNuclData = 
    new G4CrossSectionElastic( new G4ComponentGGNuclNuclXsc() );

  G4LMsdGenerator* diffGen = 0;
  if(fDiffraction) {
    diffGen = new G4LMsdGenerator("LMsdDiffraction");
    diffRatio = new G4DiffElasticRatio();
    G4AutoDelete::Register(diffRatio);    
  }

  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();
  while( (*myParticleIterator)() ) {

    G4ParticleDefinition* particle = myParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String pname = particle->GetParticleName();

    if ( pname == "anti_lambda"  || 
         pname == "anti_sigma-"  ||
         pname == "anti_sigma0"  || 
         pname == "anti_sigma+"  || 
         pname == "anti_xi-"     || 
         pname == "anti_xi0"     ||
         pname == "anti_omega-"
       ) {
      G4HadronElasticProcess* hel = new G4HadronElasticProcess();
      hel->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet( G4ChipsAntiBaryonElasticXS::Default_Name() ) );
      hel->RegisterMe( chips1 );
      pmanager->AddDiscreteProcess( hel );
      if ( verbose > 1 ) {
	G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }
      
    } else if ( pname == "lambda"  || 
                pname == "sigma-"  ||
                pname == "sigma0"  || 
                pname == "sigma+"  || 
                pname == "xi-"     || 
                pname == "xi0"     ||
                pname == "omega-"
              ) {
      G4HadronElasticProcess* hel = new G4HadronElasticProcess();
      hel->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet( G4ChipsHyperonElasticXS::Default_Name() ) );
      hel->RegisterMe( chips1 );
      pmanager->AddDiscreteProcess( hel );
      if ( verbose > 1 ) {
	G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if ( pname == "proton" ) {   
      G4HadronElasticProcess* hel = new G4HadronElasticProcess();
      hel->AddDataSet( new G4BGGNucleonElasticXS( particle ) );
      // To preserve reproducibility, a different instance of
      // G4DiffuseElastic must be used for each particle type.
      G4DiffuseElastic* protonDiffuseElastic = new G4DiffuseElastic();
      protonDiffuseElastic->SetMinEnergy( elimitDiffuse );
      hel->RegisterMe( chipsH );  // Use Chips only for Hydrogen element
      hel->RegisterMe( protonDiffuseElastic );
      pmanager->AddDiscreteProcess( hel );
      if(fDiffraction) { hel->SetDiffraction(diffGen, diffRatio); }
      if ( verbose > 1 ) {
	G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if ( pname == "neutron" ) {   
      G4HadronElasticProcess* hel = new G4HadronElasticProcess();
      hel->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronElasticXS::Default_Name()) );
      // To preserve reproducibility, a different instance of
      // G4DiffuseElastic must be used for each particle type.
      G4DiffuseElastic* neutronDiffuseElastic = new G4DiffuseElastic();
      neutronDiffuseElastic->SetMinEnergy( elimitDiffuse );
      hel->RegisterMe( chipsH );  // Use Chips only for Hydrogen element
      hel->RegisterMe( neutronDiffuseElastic );
      pmanager->AddDiscreteProcess( hel );
      if(fDiffraction) { hel->SetDiffraction(diffGen, diffRatio); }
      if ( verbose > 1 ) {
	G4cout << "### HadronElasticPhysics: " 
	       << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if ( pname == "pi-" ) { 
      G4HadronElasticProcess* hel = new G4HadronElasticProcess();
      hel->AddDataSet( new G4BGGPionElasticXS( particle ) );
      // To preserve reproducibility, a different instance of
      // G4DiffuseElastic must be used for each particle type.
      G4DiffuseElastic* pionMinusDiffuseElastic = new G4DiffuseElastic();
      pionMinusDiffuseElastic->SetMinEnergy( elimitDiffuse );
      hel->RegisterMe( chipsH );  // Use Chips only for Hydrogen element
      hel->RegisterMe( pionMinusDiffuseElastic );
      pmanager->AddDiscreteProcess( hel );
      if(fDiffraction) { hel->SetDiffraction(diffGen, diffRatio); }
      if ( verbose > 1 ) {
	G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if ( pname == "pi+" ) { 
      G4HadronElasticProcess* hel = new G4HadronElasticProcess();
      hel->AddDataSet( new G4BGGPionElasticXS( particle ) );
      // To preserve reproducibility, a different instance of
      // G4DiffuseElastic must be used for each particle type.
      G4DiffuseElastic* pionPlusDiffuseElastic = new G4DiffuseElastic();
      hel->RegisterMe( chipsH );  // Use Chips only for Hydrogen element
      pionPlusDiffuseElastic->SetMinEnergy( elimitDiffuse );
      hel->RegisterMe( pionPlusDiffuseElastic );
      pmanager->AddDiscreteProcess( hel );
      if(fDiffraction) { hel->SetDiffraction(diffGen, diffRatio); }
      if ( verbose > 1 ) {
	G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if ( pname == "kaon-"     || 
	        pname == "kaon+"     || 
	        pname == "kaon0S"    || 
	        pname == "kaon0L" 
	      ) {
      G4HadronElasticProcess* hel = new G4HadronElasticProcess();
      //AR-14Aug2017 : Replaced Chips elastic kaon cross sections with
      //               Grichine's Glauber-Gribov ones. In this way, the
      //               total (elastic + inelastic) kaon cross sections
      //               are consistent with the PDG ones.
      //               For the time being, kept Chips elastic as
      //               final-state model.
      hel->AddDataSet( theComponentGGHadronNucleusData );
      hel->RegisterMe( chips1 );
      pmanager->AddDiscreteProcess( hel );
      if(fDiffraction) { hel->SetDiffraction(diffGen, diffRatio); }
      if ( verbose > 1 ) {
	G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if ( 
                pname == "deuteron"  ||
                pname == "triton"    ||
                pname == "He3"       ||
                pname == "alpha"
              ) {
      G4HadronElasticProcess* hel = new G4HadronElasticProcess();
      hel->AddDataSet( theComponentGGNuclNuclData );
      hel->RegisterMe( diffuseNuclNuclElastic );
      pmanager->AddDiscreteProcess( hel );
      if ( verbose > 1 ) {
	G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if ( pname == "anti_proton"  ||  pname == "anti_neutron" ) {
      G4HadronElasticProcess* hel = new G4HadronElasticProcess();
      hel->AddDataSet( anucxs );
      hel->RegisterMe( chips2 );
      hel->RegisterMe( anuc );
      pmanager->AddDiscreteProcess( hel );
      if ( verbose > 1 ) {
	G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if ( pname == "anti_deuteron"  ||
                pname == "anti_triton"    ||
                pname == "anti_He3"       ||
                pname == "anti_alpha"
              ) {
      G4HadronElasticProcess* hel = new G4HadronElasticProcess();
      hel->AddDataSet( anucxs );
      hel->RegisterMe( lhep );
      hel->RegisterMe( anuc );
      pmanager->AddDiscreteProcess( hel );
      if ( verbose > 1 ) {
	G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
	       << " added for " << particle->GetParticleName() << G4endl;
      }

    } else if ( pname == "GenericIon" ) {
      G4HadronElasticProcess* hel = new G4HadronElasticProcess();
      hel->AddDataSet( theComponentGGNuclNuclData );
      hel->RegisterMe( diffuseNuclNuclElastic );
      pmanager->AddDiscreteProcess( hel );
      if ( verbose > 1 ) {
        G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
               << " added for " << particle->GetParticleName() << G4endl;
      }

    }

  }

}
