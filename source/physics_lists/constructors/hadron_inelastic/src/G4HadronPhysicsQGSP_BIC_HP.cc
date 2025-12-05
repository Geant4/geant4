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
// ClassName:   G4HadronPhysicsQGSP_BIC_HP
//
// Author: 2006 G.Folger
//
// Based on G4HadronPhysicsQGSP_BIC
//
// Modified:
// 25.04.2007 G.Folger: Add code for quasielastic
// 31.10.2012 A.Ribon: Use G4MiscBuilder
// 19.03.2013 A.Ribon: Replace LEP with FTFP and BERT
// 05.05.2020 A.Ribon: Use QGSP for antibaryons at high energies
// 07.05.2020 A.Ribon: Use QGSP for hyperons (and anti-hyperons) at high energies
// 20.05.2020 A.Ribon: Refactoring of the class (keeping same functionalities)
//
//----------------------------------------------------------------------------

#include <iomanip>   
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4NeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4QGSPNeutronBuilder.hh"
#include "G4BinaryNeutronBuilder.hh"
#include "G4BertiniNeutronBuilder.hh"
#include "G4NeutronRadCaptureHP.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4ParticleHPInelastic.hh"
#include "G4ParticleHPInelasticData.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4NeutronHPFissionData.hh"
#include "G4NeutronHPFission.hh"
#include "G4ProcessManager.hh"
#include "G4PhysListUtil.hh"
#include "G4HadronicParameters.hh"
#include "G4NuDEXNeutronCaptureModel.hh"
// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY( G4HadronPhysicsQGSP_BIC_HP );


G4HadronPhysicsQGSP_BIC_HP::G4HadronPhysicsQGSP_BIC_HP(G4int verb)
  : G4HadronPhysicsQGSP_BIC_HP( "hInelastic QGSP_BIC_HP", true )
{
  G4HadronicParameters::Instance()->SetVerboseLevel(verb);
}

G4HadronPhysicsQGSP_BIC_HP::G4HadronPhysicsQGSP_BIC_HP( const G4String& name, G4bool quasiElastic )
  :  G4HadronPhysicsQGSP_BIC( name, quasiElastic )
{
  minBIC_neutron = 19.9*CLHEP::MeV;
  G4HadronicParameters::Instance()->SetUseRFilesForXS(false);
}

void G4HadronPhysicsQGSP_BIC_HP::Neutron() {
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  G4ParticleDefinition* neutron = G4Neutron::Neutron();
  auto inel = new G4HadronInelasticProcess( "neutronInelastic", neutron );
  neutron->GetProcessManager()->AddDiscreteProcess( inel );

  G4QGSPNeutronBuilder qgs( QuasiElasticQGS );
  qgs.SetMinEnergy( minQGSP_neutron );
  qgs.Build( inel );

  G4FTFPNeutronBuilder ftf( QuasiElasticFTF );
  ftf.SetMinEnergy(minFTFP_neutron);
  ftf.SetMaxEnergy(maxFTFP_neutron);
  ftf.Build( inel );

  if ( maxBERT_neutron >  minBERT_neutron) { 
    G4BertiniNeutronBuilder bert;
    bert.SetMinEnergy( minBERT_neutron );
    bert.SetMaxEnergy( maxBERT_neutron );
    bert.Build( inel );
  }

  if ( maxBIC_neutron > 0.0 ) {
    G4BinaryNeutronBuilder bic;
    bic.SetMinEnergy( minBIC_neutron );
    bic.SetMaxEnergy( maxBIC_neutron );
    bic.Build( inel );
  }

  auto xsinel = new G4NeutronInelasticXS();
  inel->AddDataSet( xsinel );

  inel->AddDataSet( new G4ParticleHPInelasticData( neutron ) );
  auto mod = new G4ParticleHPInelastic( neutron, "NeutronHPInelastic" );
  mod->SetMaxEnergy( 20*CLHEP::MeV );
  inel->RegisterMe( mod );

  if ( useFactorXS )
    inel->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );
 
  auto capture = new G4NeutronCaptureProcess( "nCaptureHP" );
  neutron->GetProcessManager()->AddDiscreteProcess(capture);
  capture->AddDataSet( new G4NeutronHPCaptureData() );
  if (param->EnableNUDEX()) {
    capture->RegisterMe( new G4NuDEXNeutronCaptureModel() );
  } else {
    capture->RegisterMe( new G4NeutronRadCaptureHP() );
  }
  
  auto fission = new G4NeutronFissionProcess( "nFissionHP" );
  neutron->GetProcessManager()->AddDiscreteProcess(fission);
  fission->RegisterMe( new G4NeutronHPFission() );
  fission->AddDataSet( new G4ParticleHPFissionData() );
}
