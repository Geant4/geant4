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
//
//---------------------------------------------------------------------------
//
// ClassName:  HadronPhysicsNuBeam 
//
// Author: Julia Yarba, FNAL/CD (2013)
//   created from (molded after) HadronPhysicsFTFP_BERT
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsNuBeam.hh"
#include "G4QGSPLundStrFragmProtonBuilder.hh"
#include "G4ProtonBuilder.hh"
#include "G4FTFPProtonBuilder.hh"
#include "G4BertiniProtonBuilder.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "G4PhysListUtil.hh"
#include "G4HadronicParameters.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsNuBeam);

G4HadronPhysicsNuBeam::G4HadronPhysicsNuBeam(G4int) :
    G4HadronPhysicsNuBeam("hInelasticNuBeam",false)
{}

G4HadronPhysicsNuBeam::G4HadronPhysicsNuBeam(const G4String& name, G4bool quasiElastic)
    :  G4HadronPhysicsFTFP_BERT(name,quasiElastic)
{
  // specific transition energies should be defined here

  //minFTFP_neutron = 4.0*GeV;
  //maxBERT_neutron = 5.0*GeV;
  minFTFP_proton = 3.0*GeV;
  maxFTFP_proton = 101.0*GeV;
  //maxBERT_proton = 3.5*GeV;
  //minFTFP_pion = minFTFP_kaon = 3.0*GeV;
  //maxBERT_pion = maxBERT_kaon = 3.5*GeV;
}

void G4HadronPhysicsNuBeam::Proton()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();

  auto pro = new G4ProtonBuilder;
  AddBuilder(pro);
  // this is the new "custom" proton builder, tentatively for NuBeam
  //
  // no need to set the min energy because it's set in the ProBuilder
  // ... and theMax will be set via Build()
  //
  // also explicitly set quasi-elastic key ON for QGS
  // (it should be OFF for FTF, controlled by QuasiElastic)
  //
  auto qgsppro = new G4QGSPLundStrFragmProtonBuilder( true );
  AddBuilder(qgsppro);
  pro->RegisterMe(qgsppro);
  //
  // standard FTFP builder, but energy range is adjusted
  //
  auto ftfppro = new G4FTFPProtonBuilder(QuasiElastic);
  AddBuilder(ftfppro);
  pro->RegisterMe(ftfppro);
  ftfppro->SetMinEnergy(minFTFP_proton);
  ftfppro->SetMaxEnergy(maxFTFP_proton);
  //
  // standard Bertini builder
  //
  auto bertpro = new G4BertiniProtonBuilder;
  AddBuilder(bertpro);
  pro->RegisterMe(bertpro);
  bertpro->SetMaxEnergy(maxBERT_proton);
  pro->Build();

  const G4ParticleDefinition* proton = G4Proton::Proton();
  G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess(proton);
  if(inel) { 
    if( useFactorXS ) inel->MultiplyCrossSectionBy( param->XSFactorNucleonInelastic() );
  }
}

void G4HadronPhysicsNuBeam::ConstructProcess()
{
  if(G4Threading::IsMasterThread()) {
      DumpBanner();
  }
  CreateModels();
}
