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

#include <iomanip>   
#include "G4HadronPhysicsQGSP_BIC_AllHP.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProtonBuilder.hh"
#include "G4FTFPProtonBuilder.hh"
#include "G4QGSPProtonBuilder.hh"
#include "G4BinaryProtonBuilder.hh"
#include "G4ProtonPHPBuilder.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessManager.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4PhysListUtil.hh"
#include "G4HadronicParameters.hh"
// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY( G4HadronPhysicsQGSP_BIC_AllHP );


G4HadronPhysicsQGSP_BIC_AllHP::G4HadronPhysicsQGSP_BIC_AllHP(G4int verb)
  :  G4HadronPhysicsQGSP_BIC_AllHP( "hInelastic QGSP_BIC_AllHP" )
{
  G4HadronicParameters::Instance()->SetVerboseLevel(verb);
}

G4HadronPhysicsQGSP_BIC_AllHP::G4HadronPhysicsQGSP_BIC_AllHP( const G4String& name, G4bool quasiElastic )
  :  G4HadronPhysicsQGSP_BIC_HP( name, quasiElastic )
{
  minBIC_proton = 190.0*CLHEP::MeV;
  maxHP_proton  = 200.0*CLHEP::MeV;
}


void G4HadronPhysicsQGSP_BIC_AllHP::Proton() {
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  auto pro = new G4ProtonBuilder;
  AddBuilder( pro );
  auto qgs = new G4QGSPProtonBuilder( QuasiElasticQGS );
  AddBuilder( qgs );
  qgs->SetMinEnergy( minQGSP_proton );
  pro->RegisterMe( qgs );
  auto ftf = new G4FTFPProtonBuilder( QuasiElasticFTF );
  AddBuilder( ftf );
  ftf->SetMinEnergy( minFTFP_proton );
  ftf->SetMaxEnergy( maxFTFP_proton );
  pro->RegisterMe( ftf );
  auto bic = new G4BinaryProtonBuilder;
  AddBuilder( bic );
  bic->SetMinEnergy( minBIC_proton );
  bic->SetMaxEnergy( maxBIC_proton );
  pro->RegisterMe( bic );
  auto hp = new G4ProtonPHPBuilder;
  AddBuilder( hp );
  hp->SetMaxEnergy( maxHP_proton );
  pro->RegisterMe( hp );
  pro->Build();

  const G4ParticleDefinition* proton = G4Proton::Proton();
  G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess(proton);
  if(nullptr != inel) { 
    if( useFactorXS ) inel->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );
  }
}
