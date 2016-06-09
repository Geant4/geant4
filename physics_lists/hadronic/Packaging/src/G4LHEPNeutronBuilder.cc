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
 #include "G4LHEPNeutronBuilder.hh"
 #include "G4ParticleDefinition.hh"
 #include "G4ParticleTable.hh"
 #include "G4ProcessManager.hh"

 G4LHEPNeutronBuilder::
 G4LHEPNeutronBuilder() 
 {
   theMin = 0;
   theIMin = 0;
 }

 G4LHEPNeutronBuilder::
 ~G4LHEPNeutronBuilder() 
{
  delete theLENeutronModel;
  delete theHENeutronModel;
}

 void G4LHEPNeutronBuilder::
 Build(G4NeutronInelasticProcess * aP)
 {
   theLENeutronModel = new G4LENeutronInelastic();
   theHENeutronModel = new G4HENeutronInelastic();
   theLENeutronModel->SetMinEnergy(theIMin);
   theLENeutronModel->SetMaxEnergy(55*GeV);
   theHENeutronModel->SetMinEnergy(25*GeV);
   aP->RegisterMe(theLENeutronModel);
   aP->RegisterMe(theHENeutronModel);
 }

 void G4LHEPNeutronBuilder::
 Build(G4HadronFissionProcess * aP)
 {
   theNeutronFissionModel = new G4LFission();
   theNeutronFissionModel->SetMinEnergy(theMin);
   aP->RegisterMe(theNeutronFissionModel);
 }

 void G4LHEPNeutronBuilder::
 Build(G4HadronElasticProcess *)
 {
 }

 void G4LHEPNeutronBuilder::
 Build(G4HadronCaptureProcess * aP)
 {
   theNeutronCaptureModel = new G4LCapture();
   theNeutronCaptureModel->SetMinEnergy(theMin);
   aP->RegisterMe(theNeutronCaptureModel);
 }

 // 2002 by J.P. Wellisch
