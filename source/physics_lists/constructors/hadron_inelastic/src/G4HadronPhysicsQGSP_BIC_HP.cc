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
#include "G4NeutronPHPBuilder.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4LFission.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessManager.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4PhysListUtil.hh"
#include "G4HadronicParameters.hh"
// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY( G4HadronPhysicsQGSP_BIC_HP );


G4HadronPhysicsQGSP_BIC_HP::G4HadronPhysicsQGSP_BIC_HP(G4int verb)
  : G4HadronPhysicsQGSP_BIC_HP( "hInelastic QGSP_BIC_HP" )
{
  G4HadronicParameters::Instance()->SetVerboseLevel(verb);
}

G4HadronPhysicsQGSP_BIC_HP::G4HadronPhysicsQGSP_BIC_HP( const G4String& name, G4bool quasiElastic )
  :  G4HadronPhysicsQGSP_BIC( name, quasiElastic )
{
  minBIC_neutron = 19.9*MeV;
}

void G4HadronPhysicsQGSP_BIC_HP::Neutron() {
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  auto neu = new G4NeutronBuilder( true );  // Fission on
  AddBuilder( neu );
  auto qgs = new G4QGSPNeutronBuilder( QuasiElasticQGS );
  AddBuilder( qgs );
  qgs->SetMinEnergy( minQGSP_neutron );
  neu->RegisterMe( qgs );
  auto ftf = new G4FTFPNeutronBuilder( QuasiElasticFTF );
  AddBuilder( ftf );
  ftf->SetMinEnergy( minFTFP_neutron );
  ftf->SetMaxEnergy( maxFTFP_neutron );
  neu->RegisterMe( ftf );
  auto bic = new G4BinaryNeutronBuilder;
  AddBuilder( bic );
  bic->SetMinEnergy( minBIC_neutron );
  bic->SetMaxEnergy( maxBIC_neutron );
  neu->RegisterMe( bic );
  auto hp = new G4NeutronPHPBuilder;
  AddBuilder( hp );
  neu->RegisterMe( hp );
  neu->Build();

  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess( neutron );
  if(inel) { 
    if( useFactorXS ) inel->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );
  }
  G4HadronicProcess* capture = G4PhysListUtil::FindCaptureProcess( neutron );
  if ( capture ) {
    G4NeutronRadCapture* theNeutronRadCapture = new G4NeutronRadCapture;
    theNeutronRadCapture->SetMinEnergy( minBIC_neutron );
    capture->RegisterMe( theNeutronRadCapture );
  }
  G4HadronicProcess* fission = G4PhysListUtil::FindFissionProcess( neutron );
  if ( fission ) {
    G4LFission* theNeutronLEPFission = new G4LFission;
    theNeutronLEPFission->SetMinEnergy( minBIC_neutron );
    theNeutronLEPFission->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
    fission->RegisterMe( theNeutronLEPFission );
  }
}
