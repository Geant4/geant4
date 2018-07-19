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
// $Id: G4FTFPAntiBarionBuilder.cc 103555 2017-04-18 09:04:37Z gcosmo $
// GEANT4 tag $Name:  $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4FTFPAntiBarionBuilder
//
// Author: 2011 J. Apostolakis
//
// Modified:
//  2011.02.22  J.Apostolakis - Started from G4FTFPPiKBuilder
//----------------------------------------------------------------------------
//
#include "G4FTFPAntiBarionBuilder.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"  // For anti-ions
#include "G4CrossSectionInelastic.hh"

// #include "G4AntiBarionNuclearCrossSection.hh"

G4FTFPAntiBarionBuilder::
G4FTFPAntiBarionBuilder(G4bool quasiElastic) 
{
  theAntiNucleonData = 
    new G4CrossSectionInelastic(theAntiNucleonXS=new G4ComponentAntiNuclNuclearXS());

  theMin =   0.0*GeV;
  theMax = 100.0*TeV;
  theModel = new G4TheoFSGenerator("FTFP");

  theStringModel = new G4FTFModel;
  theStringDecay = new G4ExcitedStringDecay(theLund = new G4LundStringFragmentation);
  theStringModel->SetFragmentationModel(theStringDecay);

  theCascade = new G4GeneratorPrecompoundInterface();

  theModel->SetHighEnergyGenerator(theStringModel);
  if (quasiElastic)
  {
     theQuasiElastic=new G4QuasiElasticChannel;
     theModel->SetQuasiElasticChannel(theQuasiElastic);
  } else 
  {  theQuasiElastic=0;}  

  theModel->SetTransport(theCascade);
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(100*TeV);
}

G4FTFPAntiBarionBuilder::~G4FTFPAntiBarionBuilder() 
{
  delete theStringDecay;
  delete theStringModel;
  //delete theModel;
  if ( theQuasiElastic ) delete theQuasiElastic;
  delete theLund;
  delete theAntiNucleonXS;
  //delete theAntiNucleonData;
}

void G4FTFPAntiBarionBuilder::
Build(G4AntiProtonInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  aP->AddDataSet(theAntiNucleonData);
  aP->RegisterMe(theModel);
}

void G4FTFPAntiBarionBuilder::
Build(G4AntiNeutronInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  aP->AddDataSet(theAntiNucleonData);
  aP->RegisterMe(theModel);
}

void G4FTFPAntiBarionBuilder::
Build(G4AntiDeuteronInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  aP->AddDataSet(theAntiNucleonData);
  aP->RegisterMe(theModel);
}

void G4FTFPAntiBarionBuilder::
Build(G4AntiTritonInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  aP->AddDataSet(theAntiNucleonData);
  aP->RegisterMe(theModel);
}

void G4FTFPAntiBarionBuilder::
Build(G4AntiHe3InelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  aP->AddDataSet(theAntiNucleonData);
  aP->RegisterMe(theModel);
}

void G4FTFPAntiBarionBuilder::
Build(G4AntiAlphaInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  aP->AddDataSet(theAntiNucleonData);
  aP->RegisterMe(theModel);
}
