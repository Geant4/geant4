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
// $Id: G4QGSPPiKBuilder.cc 103555 2017-04-18 09:04:37Z gcosmo $
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
#include "G4PiNuclearCrossSection.hh"
#include "G4CrossSectionPairGG.hh"
#include "G4ChipsKaonMinusInelasticXS.hh"
#include "G4ChipsKaonPlusInelasticXS.hh"
#include "G4ChipsKaonZeroInelasticXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"


G4QGSPPiKBuilder::
G4QGSPPiKBuilder(G4bool quasiElastic) 
{
  thePiData = new G4CrossSectionPairGG((G4PiNuclearCrossSection*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4PiNuclearCrossSection::Default_Name()), 91*GeV);
    
  theMin = 12*GeV;
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
  theModel->SetMaxEnergy(100*TeV);
  aP->AddDataSet(thePiData);
  aP->RegisterMe(theModel);
}

void G4QGSPPiKBuilder::
Build(G4PionMinusInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(100*TeV);
  aP->AddDataSet(thePiData);
  aP->RegisterMe(theModel);
}

void G4QGSPPiKBuilder::
Build(G4KaonPlusInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(100*TeV);
    aP->RegisterMe(theModel);
    aP->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonPlusInelasticXS::Default_Name()));
}

void G4QGSPPiKBuilder::
Build(G4KaonMinusInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(100*TeV);
  aP->RegisterMe(theModel);
  aP->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonMinusInelasticXS::Default_Name()));
}

void G4QGSPPiKBuilder::
Build(G4KaonZeroLInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(100*TeV);
  aP->RegisterMe(theModel);
  aP->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
}

void G4QGSPPiKBuilder::
Build(G4KaonZeroSInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(100*TeV);
  aP->RegisterMe(theModel);
}

// 2002 by J.P. Wellisch
