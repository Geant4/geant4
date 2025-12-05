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
// ClassName:   G4HadronPhysicsQGS_BIC
//
// Author: 2007 Gunter Folger
//     created from G4HadronPhysicsQGSP_BIC  by H.P.Wellisch
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4HadronPhysicsQGS_BIC.hh"
#include "G4PionBuilder.hh"
#include "G4BinaryPionBuilder.hh"
#include "G4BertiniPionBuilder.hh"
#include "G4FTFBinaryPionBuilder.hh"
#include "G4QGSBinaryPionBuilder.hh"

#include "G4KaonBuilder.hh"
#include "G4BertiniKaonBuilder.hh"
#include "G4FTFBinaryKaonBuilder.hh"
#include "G4QGSBinaryKaonBuilder.hh"

#include "G4BertiniProtonBuilder.hh"
#include "G4FTFBinaryProtonBuilder.hh"
#include "G4QGSBinaryProtonBuilder.hh"
#include "G4BinaryProtonBuilder.hh"

#include "G4BertiniNeutronBuilder.hh"
#include "G4FTFBinaryNeutronBuilder.hh"
#include "G4QGSBinaryNeutronBuilder.hh"
#include "G4BinaryNeutronBuilder.hh"
#include "G4ParticleInelasticXS.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4NeutronRadCapture.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"

#include "G4PhysListUtil.hh"
#include "G4HadParticles.hh"
#include "G4HadronicParameters.hh"
#include "G4ProcessManager.hh"

#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGS_BIC);

G4HadronPhysicsQGS_BIC::G4HadronPhysicsQGS_BIC(G4int verb)
    : G4HadronPhysicsQGS_BIC("hInelastic QGS_BIC",true) 
{
  G4HadronicParameters::Instance()->SetVerboseLevel(verb);
}

G4HadronPhysicsQGS_BIC::G4HadronPhysicsQGS_BIC(const G4String& name, G4bool qe)
  : G4HadronPhysicsQGSP_BERT(name, qe) 
{
  maxBIC_proton = maxBIC_neutron  = 1.5*CLHEP::GeV;
  minBERT_proton = minBERT_neutron  = 1.0*CLHEP::GeV;
}

void G4HadronPhysicsQGS_BIC::Neutron()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  auto inel = new G4HadronInelasticProcess( "neutronInelastic", neutron );
  neutron->GetProcessManager()->AddDiscreteProcess( inel );

  G4QGSBinaryNeutronBuilder qgs( QuasiElasticQGS );
  qgs.SetMinEnergy( minQGSP_neutron );
  qgs.Build( inel );

  G4FTFBinaryNeutronBuilder ftf( QuasiElasticFTF );
  ftf.SetMinEnergy( minFTFP_neutron );
  ftf.SetMaxEnergy( maxFTFP_neutron );
  ftf.Build( inel );
  
  G4BertiniNeutronBuilder bert;
  bert.SetMinEnergy( minBERT_neutron );
  bert.SetMaxEnergy( maxBERT_neutron );
  bert.Build( inel );

  if ( maxBIC_neutron > 0.0 ) {
    G4BinaryNeutronBuilder bic;
    bic.SetMaxEnergy( maxBIC_neutron );
    bic.Build( inel );
  }

  inel->AddDataSet( new G4NeutronInelasticXS() ); 
  if ( useFactorXS ) {
    inel->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );
  }
  auto capture = new G4NeutronCaptureProcess( "nCaptureXS" );
  neutron->GetProcessManager()->AddDiscreteProcess(capture);
  capture->AddDataSet( new G4NeutronCaptureXS() );
  capture->RegisterMe( new G4NeutronRadCapture() );
}

void G4HadronPhysicsQGS_BIC::Proton()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  const G4ParticleDefinition* proton = G4Proton::Proton();
  auto inel = new G4HadronInelasticProcess( "protonInelastic", proton );
  proton->GetProcessManager()->AddDiscreteProcess( inel );

  G4QGSBinaryProtonBuilder qgs(QuasiElasticQGS);
  qgs.SetMinEnergy(minQGSP_proton);
  qgs.Build( inel );

  G4FTFBinaryProtonBuilder ftf(QuasiElasticFTF);
  ftf.SetMinEnergy( minFTFP_proton );
  ftf.SetMaxEnergy( maxFTFP_proton );
  ftf.Build( inel );

  G4BertiniProtonBuilder bert;
  bert.SetMinEnergy( minBERT_proton );
  bert.SetMaxEnergy( maxBERT_proton );
  bert.Build( inel );

  if ( maxBIC_proton > 0.0 ) {
    G4BinaryProtonBuilder bic;
    bic.SetMaxEnergy( maxBIC_proton);
    bic.Build( inel );
  }

  auto xsinel = new G4ParticleInelasticXS( proton );
  inel->AddDataSet( xsinel );

  if ( useFactorXS ) {
    inel->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );
  }
}
