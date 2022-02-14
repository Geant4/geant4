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
//--------------------------------------------------------------------------
// ClassName: G4QGSPAntiBarionBuilder
// Author: Alberto Ribon
// Date: May 2020
//
// Modified:
//---------------------------------------------------------------------------

#include "G4QGSPAntiBarionBuilder.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4HadronicParameters.hh"
#include "G4TheoFSGenerator.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4QuasiElasticChannel.hh"


G4QGSPAntiBarionBuilder::G4QGSPAntiBarionBuilder( G4bool quasiElastic ) {
  G4CrossSectionDataSetRegistry* xsreg = G4CrossSectionDataSetRegistry::Instance();
  G4VComponentCrossSection* theAntiNucleonXS = xsreg->GetComponentCrossSection( "AntiAGlauber" );
  if ( ! theAntiNucleonXS ) theAntiNucleonXS = new G4ComponentAntiNuclNuclearXS;
  theAntiNucleonData = new G4CrossSectionInelastic( theAntiNucleonXS );
  theMin = G4HadronicParameters::Instance()->GetMinEnergyTransitionQGS_FTF();  // Safe choice
  theMax = G4HadronicParameters::Instance()->GetMaxEnergy();
  // The main model, QGSP, applicable only for anti_proton and anti_neutron  
  theQGSmodel = new G4TheoFSGenerator( "QGSP" );
  G4QGSModel< G4QGSParticipants >* theStringModel = new G4QGSModel< G4QGSParticipants >;
  G4ExcitedStringDecay* theStringDecay = new G4ExcitedStringDecay( new G4QGSMFragmentation );
  theStringModel->SetFragmentationModel( theStringDecay );
  G4GeneratorPrecompoundInterface* theCascade = new G4GeneratorPrecompoundInterface();
  theQGSmodel->SetTransport( theCascade );
  theQGSmodel->SetHighEnergyGenerator( theStringModel );
  if ( quasiElastic ) theQGSmodel->SetQuasiElasticChannel( new G4QuasiElasticChannel );
  theQGSmodel->SetTransport(theCascade);
  theQGSmodel->SetMinEnergy( theMin );
  theQGSmodel->SetMaxEnergy( theMax );
  // The auxilary model, FTFP, needed for anti_deuteron, anti_triton, anti_He3 and anti_alpha
  theFTFmodel = new G4TheoFSGenerator( "FTFP" );
  G4FTFModel* theStringModel2 = new G4FTFModel;
  theStringModel2->SetFragmentationModel( new G4ExcitedStringDecay );
  G4GeneratorPrecompoundInterface* theCascade2 = new G4GeneratorPrecompoundInterface;
  theFTFmodel->SetHighEnergyGenerator( theStringModel2 );
  G4double quasiElasticFTF = false;  // Use built-in quasi-elastic (not add-on)
  if ( quasiElasticFTF ) theFTFmodel->SetQuasiElasticChannel( new G4QuasiElasticChannel );
  theFTFmodel->SetTransport( theCascade2 );
  theFTFmodel->SetMinEnergy( theMin );
  theFTFmodel->SetMaxEnergy( theMax );
}


void G4QGSPAntiBarionBuilder::Build( G4HadronInelasticProcess* aP ) {
  if ( aP->GetParticleDefinition()  &&  aP->GetParticleDefinition()->GetBaryonNumber() < -1 ) {
    // Light anti-ions: for the time being QGSP cannot be applied, use FTFP
    theFTFmodel->SetMinEnergy( theMin );
    theFTFmodel->SetMaxEnergy( theMax );
    aP->RegisterMe( theFTFmodel );
  } else {
    // Anti-proton and anti-neutron: use QGSP
    theQGSmodel->SetMinEnergy( theMin );
    theQGSmodel->SetMaxEnergy( theMax );
    aP->RegisterMe( theQGSmodel );
  }
  aP->AddDataSet( theAntiNucleonData );
}
