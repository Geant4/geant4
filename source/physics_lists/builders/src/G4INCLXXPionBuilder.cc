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
 #include "G4INCLXXPionBuilder.hh"
 #include "G4SystemOfUnits.hh"
 #include "G4ParticleDefinition.hh"
 #include "G4ParticleTable.hh"
 #include "G4ProcessManager.hh"
 #include "G4PiNuclearCrossSection.hh"

 G4INCLXXPionBuilder::
 G4INCLXXPionBuilder() 
 {
   theMin = 0*GeV;
   theMax = 20.0*GeV;
   theModel = new G4INCLXXInterface();
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(theMax); 
   thePiData = new G4PiNuclearCrossSection();
 }

 void G4INCLXXPionBuilder::
 Build(G4PionPlusInelasticProcess * aP)
 {
   aP->RegisterMe(theModel);
   aP->AddDataSet(thePiData);
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(theMax);
 }

 void G4INCLXXPionBuilder::
 Build(G4PionMinusInelasticProcess * aP)
 {
   aP->RegisterMe(theModel);
   aP->AddDataSet(thePiData);
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(theMax);
 }


 // 2011 by P. Kaitaniemi

