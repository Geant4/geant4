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
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronPhysicsQGSP_BERT_HP
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 15.12.2005 G.Folger: migration to non static particles
// 08.06.2006 V.Ivanchenko: remove stopping
// 20.06.2006 G.Folger: Bertini applies to Kaons, i.e. use SetMinEnergy
//            instead of SetMinPionEnergy
// 25.04.2007 G.Folger: Add code for quasielastic
// 31.10.2012 A.Ribon: Use G4MiscBuilder
// 19.03.2013 A.Ribon: Replace LEP with FTFP
// 12.10.2023 V.Ivanchenko added usage of alternative neutron
//            HP model and cross section
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsQGSP_BERT_HP.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4NeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4QGSPNeutronBuilder.hh"
#include "G4BertiniNeutronBuilder.hh"
#include "G4NeutronPHPBuilder.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4NeutronCaptureProcess.hh"
#include "G4NeutronFissionProcess.hh"

#include "G4NeutronRadCaptureHP.hh"
#include "G4NeutronHPCaptureXS.hh"
#include "G4NeutronHPFissionXS.hh"
#include "G4NeutronHPInelasticXS.hh"
#include "G4NeutronHPInelasticVI.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4ParticleHPInelastic.hh"
#include "G4NeutronFissionVI.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessManager.hh"

#include "G4HadronicParameters.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGSP_BERT_HP);

G4HadronPhysicsQGSP_BERT_HP::G4HadronPhysicsQGSP_BERT_HP(G4int verb)
    :  G4HadronPhysicsQGSP_BERT_HP("hInelastic QGSP_BERT_HP")
{
  G4HadronicParameters::Instance()->SetVerboseLevel(verb);
}

G4HadronPhysicsQGSP_BERT_HP::G4HadronPhysicsQGSP_BERT_HP(const G4String& name, G4bool /*quasiElastic */ )
    :  G4HadronPhysicsQGSP_BERT(name)
{
  minBERT_neutron = 19.9*MeV;
  auto param = G4HadronicParameters::Instance();
  // HP is inconsistent with the neutron general process
  param->SetEnableNeutronGeneralProcess(false);
}

void G4HadronPhysicsQGSP_BERT_HP::Neutron()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  auto inel = new G4HadronInelasticProcess( "neutronInelastic", neutron );
  neutron->GetProcessManager()->AddDiscreteProcess(inel);

  G4QGSPNeutronBuilder qgs(QuasiElasticQGS);
  qgs.SetMinEnergy(minQGSP_neutron);
  qgs.Build(inel);

  G4FTFPNeutronBuilder ftf(QuasiElasticFTF);
  ftf.SetMinEnergy(minFTFP_neutron);
  ftf.SetMaxEnergy(maxFTFP_neutron);
  ftf.Build(inel);

  G4BertiniNeutronBuilder bert;
  bert.SetMinEnergy(minBERT_neutron);
  bert.SetMaxEnergy(maxBERT_neutron);
  bert.Build(inel);

  auto xsinel = new G4NeutronInelasticXS();
  inel->AddDataSet( xsinel );
  inel->AddDataSet( new G4NeutronHPInelasticXS() );
  auto mod = new G4NeutronHPInelasticVI();
  mod->SetMaxEnergy( 20*CLHEP::MeV );
  inel->RegisterMe( mod );
  if ( useFactorXS )
    inel->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );
 
  auto capture = new G4NeutronCaptureProcess( "nCaptureHP" );
  neutron->GetProcessManager()->AddDiscreteProcess(capture);
  capture->AddDataSet( new G4NeutronHPCaptureXS() );
  capture->RegisterMe( new G4NeutronRadCaptureHP() );
  
  auto fission = new G4NeutronFissionProcess( "nFissionHP" );
  neutron->GetProcessManager()->AddDiscreteProcess(fission);
  fission->AddDataSet( new G4NeutronHPFissionXS() );
  fission->RegisterMe( new G4NeutronFissionVI() );
}
