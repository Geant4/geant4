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
/// \file hadronic/Hadr02/src/HadronPhysicsCRMC_FTFP_BERT.cc
/// \brief Implementation of the CRMC_FTFP_BERT class methods
//
//
//---------------------------------------------------------------------------
//
// ClassName:   HadronPhysicsCRMC_FTFP_BERT
//
// Authors: 2018 Alberto Ribon
//
// Modified:
// -  18-May-2021 Alberto Ribon : Migrated to newer physics constructor
//                                and used the latest Geant4-CRMC interface.
//
//----------------------------------------------------------------------------
//
#ifdef G4_USE_CRMC

#include <iomanip>   
#include "HadronPhysicsCRMC_FTFP_BERT.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PionBuilder.hh"
#include "G4KaonBuilder.hh"
#include "G4BertiniPionBuilder.hh"
#include "G4BertiniKaonBuilder.hh"
#include "G4FTFPPionBuilder.hh"
#include "G4FTFPKaonBuilder.hh"
#include "CRMCPionBuilder.hh"
#include "CRMCKaonBuilder.hh"
#include "G4ProtonBuilder.hh"
#include "G4BertiniProtonBuilder.hh"
#include "G4FTFPProtonBuilder.hh"
#include "CRMCProtonBuilder.hh"
#include "G4NeutronBuilder.hh"
#include "G4BertiniNeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "CRMCNeutronBuilder.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4PhysListUtil.hh"
#include "G4HadParticles.hh"
#include "G4PhysicsConstructorFactory.hh"

G4_DECLARE_PHYSCONSTR_FACTORY( HadronPhysicsCRMC_FTFP_BERT );

const std::array< std::string, 13 > HadronPhysicsCRMC_FTFP_BERT::fModelNames = {
  "EPOS-LHC", "EPOS-1.99", "QGSJET-01", "", "", "",
  "SIBYLL-2.3", "QGSJETII-04", "", "", "", "QGSJETII-03", "DPMJET-3.06" };          


HadronPhysicsCRMC_FTFP_BERT::HadronPhysicsCRMC_FTFP_BERT( G4int )
  : HadronPhysicsCRMC_FTFP_BERT( "hInelastic CRMC_FTFP_BERT", false ) {}


HadronPhysicsCRMC_FTFP_BERT::HadronPhysicsCRMC_FTFP_BERT( const G4String& name, G4bool qe )
  : G4HadronPhysicsFTFP_BERT( name, qe ) {
  fModel   = 0;          //***LOOKHERE*** CRMC model: 0:EPOS-LHC, 1:EPOS-1.99, 2:QGSJET:01, 6:SIBYLL-2.3,
                        //                           7:QGSJETII-04, 11:QGSJETII-03, 12:DPMJET-3.06
  fMinCRMC = 100.0*GeV;  //***LOOKHERE*** CRMC model is applied only above this projectile lab energy
  fMaxFTFP = 110.0*GeV;  //***LOOKHERE*** FTFP model is applied only below this projectile lab energy
}


HadronPhysicsCRMC_FTFP_BERT::~HadronPhysicsCRMC_FTFP_BERT() {} 


void HadronPhysicsCRMC_FTFP_BERT::Neutron() {
  auto neutronBuilder = new G4NeutronBuilder;
  AddBuilder( neutronBuilder );
  auto ftfpnBuilder = new G4FTFPNeutronBuilder( QuasiElastic );
  ftfpnBuilder->SetMinEnergy( minFTFP_neutron );
  ftfpnBuilder->SetMaxEnergy( fMaxFTFP );
  AddBuilder( ftfpnBuilder );
  neutronBuilder->RegisterMe( ftfpnBuilder );
  auto bertnBuilder = new G4BertiniNeutronBuilder;
  bertnBuilder->SetMaxEnergy( maxBERT_neutron );
  AddBuilder( bertnBuilder );
  neutronBuilder->RegisterMe( bertnBuilder );
  auto crmcnBuilder = new CRMCNeutronBuilder( fModel, fModelNames[fModel] );
  crmcnBuilder->SetMinEnergy( fMinCRMC );
  AddBuilder( crmcnBuilder );
  neutronBuilder->RegisterMe( crmcnBuilder );
  neutronBuilder->Build();
  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  G4HadronicProcess* inel = G4PhysListUtil::FindInelasticProcess( neutron );
  if ( inel ) inel->AddDataSet( new G4NeutronInelasticXS ); 
  G4HadronicProcess* capture = G4PhysListUtil::FindCaptureProcess( neutron );
  if ( capture ) capture->RegisterMe( new G4NeutronRadCapture );
}


void HadronPhysicsCRMC_FTFP_BERT::Proton() {
  auto protonBuilder = new G4ProtonBuilder;
  AddBuilder( protonBuilder );
  auto ftfppBuilder = new G4FTFPProtonBuilder( QuasiElastic );
  ftfppBuilder->SetMinEnergy( minFTFP_proton );
  ftfppBuilder->SetMaxEnergy( fMaxFTFP );
  AddBuilder( ftfppBuilder );
  protonBuilder->RegisterMe( ftfppBuilder );
  auto bertpBuilder = new G4BertiniProtonBuilder;
  bertpBuilder->SetMaxEnergy( maxBERT_proton );
  AddBuilder( bertpBuilder );
  protonBuilder->RegisterMe( bertpBuilder );
  auto crmcpBuilder = new CRMCProtonBuilder( fModel, fModelNames[fModel] );
  crmcpBuilder->SetMinEnergy( fMinCRMC );
  AddBuilder( crmcpBuilder );
  protonBuilder->RegisterMe( crmcpBuilder );
  protonBuilder->Build();
} 


void HadronPhysicsCRMC_FTFP_BERT::Pion() {
  auto pionBuilder = new G4PionBuilder;
  AddBuilder( pionBuilder );
  auto ftfppiBuilder = new G4FTFPPionBuilder( QuasiElastic );
  ftfppiBuilder->SetMinEnergy( minFTFP_pion );
  ftfppiBuilder->SetMaxEnergy( fMaxFTFP );
  AddBuilder( ftfppiBuilder );
  pionBuilder->RegisterMe( ftfppiBuilder );
  auto bertpiBuilder = new G4BertiniPionBuilder;
  bertpiBuilder->SetMaxEnergy( maxBERT_pion );
  AddBuilder( bertpiBuilder );
  pionBuilder->RegisterMe( bertpiBuilder );
  auto crmcpiBuilder = new CRMCPionBuilder( fModel, fModelNames[fModel] );
  crmcpiBuilder->SetMinEnergy( fMinCRMC );
  AddBuilder( crmcpiBuilder );
  pionBuilder->RegisterMe( crmcpiBuilder );
  pionBuilder->Build();
}


void HadronPhysicsCRMC_FTFP_BERT::Kaon() {
  auto kaonBuilder = new G4KaonBuilder;
  AddBuilder( kaonBuilder );
  auto ftfpkBuilder = new G4FTFPKaonBuilder( QuasiElastic );
  ftfpkBuilder->SetMinEnergy( minFTFP_kaon );
  ftfpkBuilder->SetMaxEnergy( fMaxFTFP );
  AddBuilder( ftfpkBuilder );
  kaonBuilder->RegisterMe( ftfpkBuilder );
  auto bertkBuilder = new G4BertiniKaonBuilder;
  bertkBuilder->SetMaxEnergy( maxBERT_kaon );
  AddBuilder( bertkBuilder );
  kaonBuilder->RegisterMe( bertkBuilder );
  auto crmckBuilder = new CRMCKaonBuilder( fModel, fModelNames[fModel] );
  crmckBuilder->SetMinEnergy( fMinCRMC );
  AddBuilder( crmckBuilder );
  kaonBuilder->RegisterMe( crmckBuilder );  
  kaonBuilder->Build();
}  

#endif //G4_USE_CRMC
