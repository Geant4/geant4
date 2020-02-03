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
// ClassName:   G4QGSPPiKBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 17.11.2010 G.Folger, use G4CrossSectionPairGG for relativistic rise of cross
//             section at high energies.
// 30.03.2009 V.Ivanchenko create cross section by new
//
//----------------------------------------------------------------------------
//
#include "G4QGSPPiKBuilder.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4BGGPionInelasticXS.hh"
#include "G4HadronicParameters.hh"


G4QGSPPiKBuilder::
G4QGSPPiKBuilder(G4bool quasiElastic) 
{    
  theMin = G4HadronicParameters::Instance()->GetMinEnergyTransitionQGS_FTF();
  theModel = new G4TheoFSGenerator("QGSP");

  theStringModel = new G4QGSModel< G4QGSParticipants >;
  theStringDecay = new G4ExcitedStringDecay(theQGSM = new G4QGSMFragmentation);
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
}
 
G4QGSPPiKBuilder::
~G4QGSPPiKBuilder() 
{
  if ( theQuasiElastic ) delete theQuasiElastic;
  delete theStringDecay;
  delete theStringModel;
  delete theQGSM;
}

void G4QGSPPiKBuilder::
Build(G4PionPlusInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  aP->AddDataSet( new G4BGGPionInelasticXS( G4PionPlus::Definition() ) );
  aP->RegisterMe(theModel);
}

void G4QGSPPiKBuilder::
Build(G4PionMinusInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  aP->AddDataSet( new G4BGGPionInelasticXS( G4PionMinus::Definition() ) );
  aP->RegisterMe(theModel);
}

void G4QGSPPiKBuilder::
Build(G4KaonPlusInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  aP->RegisterMe(theModel);
}

void G4QGSPPiKBuilder::
Build(G4KaonMinusInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  aP->RegisterMe(theModel);
}

void G4QGSPPiKBuilder::
Build(G4KaonZeroLInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  aP->RegisterMe(theModel);
}

void G4QGSPPiKBuilder::
Build(G4KaonZeroSInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  aP->RegisterMe(theModel);
}

