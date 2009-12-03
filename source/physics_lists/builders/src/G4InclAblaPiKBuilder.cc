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
 #include "G4InclAblaPiKBuilder.hh"
 #include "G4ParticleDefinition.hh"
 #include "G4ParticleTable.hh"
 #include "G4ProcessManager.hh"

 G4InclAblaPiKBuilder::
 G4InclAblaPiKBuilder() 
 {
   theMin = 0*GeV;
   theMax = 3.0*GeV;
   theModel = new G4InclAblaCascadeInterface;
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(theMax); 
   theBertiniModel = new G4CascadeInterface;
   theBertiniModel->SetMinEnergy(theMin);
   theBertiniModel->SetMaxEnergy(theMax); 
 }

 G4InclAblaPiKBuilder::
 ~G4InclAblaPiKBuilder() 
{
  delete theModel;
}

 void G4InclAblaPiKBuilder::
 Build(G4PionPlusInelasticProcess * aP)
 {
   aP->RegisterMe(theModel);
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(theMax);
 }

 void G4InclAblaPiKBuilder::
 Build(G4PionMinusInelasticProcess * aP)
 {
   aP->RegisterMe(theModel);
   aP->AddDataSet(&thePiData);
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(theMax);
 }

 void G4InclAblaPiKBuilder::
 Build(G4HadronElasticProcess * ) {}

 void G4InclAblaPiKBuilder::
 Build(G4KaonPlusInelasticProcess *aP)
 {
   aP->RegisterMe(theBertiniModel);
   theBertiniModel->SetMinEnergy(theMin);
   theBertiniModel->SetMaxEnergy(theMax);
 }

 void G4InclAblaPiKBuilder::
 Build(G4KaonMinusInelasticProcess *aP)
 {
   aP->RegisterMe(theBertiniModel);
   theBertiniModel->SetMinEnergy(theMin);
   theBertiniModel->SetMaxEnergy(theMax);
 }

 void G4InclAblaPiKBuilder::
 Build(G4KaonZeroLInelasticProcess *aP)
 {
   aP->RegisterMe(theBertiniModel);
   theBertiniModel->SetMinEnergy(theMin);
   theBertiniModel->SetMaxEnergy(theMax);
 }

 void G4InclAblaPiKBuilder::
 Build(G4KaonZeroSInelasticProcess *aP)
 {
   aP->RegisterMe(theBertiniModel);
   theBertiniModel->SetMinEnergy(theMin);
   theBertiniModel->SetMaxEnergy(theMax);
 }

 // 2002 by J.P. Wellisch
