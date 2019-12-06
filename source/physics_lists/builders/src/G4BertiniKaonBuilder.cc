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
// ClassName:   G4BertiniKaonBuilder
//
// Author: 2011 G.Folger
//  devired from G4BertiniPionBuilder
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include "G4BertiniKaonBuilder.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4CrossSectionInelastic.hh"

#include "G4HadronicParameters.hh"


G4BertiniKaonBuilder::
G4BertiniKaonBuilder() 
 {
   KaonMinusCrossSection = KaonPlusCrossSection = KaonZeroCrossSection = 
     new G4CrossSectionInelastic( new G4ComponentGGHadronNucleusXsc );
   theMin = 0.0;
   theMax = G4HadronicParameters::Instance()->GetMaxEnergyTransitionFTF_Cascade();
   theModel = new G4CascadeInterface;
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(theMax); 
 }

void G4BertiniKaonBuilder::
Build(G4KaonPlusInelasticProcess * aP)
 {
   aP->RegisterMe(theModel);
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(theMax);
   aP->AddDataSet(KaonPlusCrossSection);
 }

void G4BertiniKaonBuilder::
Build(G4KaonMinusInelasticProcess * aP)
 {
   aP->RegisterMe(theModel);
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(theMax);
   aP->AddDataSet(KaonMinusCrossSection);
 }

void G4BertiniKaonBuilder::
Build(G4KaonZeroLInelasticProcess * aP)
 {
   aP->RegisterMe(theModel);
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(theMax);
   aP->AddDataSet(KaonZeroCrossSection);
 }

void G4BertiniKaonBuilder::
Build(G4KaonZeroSInelasticProcess * aP)
 {
   aP->RegisterMe(theModel);
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(theMax);
   aP->AddDataSet(KaonZeroCrossSection);
 }

