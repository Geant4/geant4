//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
 #include "G4QGSPProtonBuilder.hh"
 #include "G4ParticleDefinition.hh"
 #include "G4ParticleTable.hh"
 #include "G4ProcessManager.hh"

 G4QGSPProtonBuilder::
 G4QGSPProtonBuilder() 
 {
   theMin = 12*GeV;
   theModel = new G4TheoFSGenerator;

   theStringModel = new G4QGSModel< G4QGSParticipants >;
   theStringDecay = new G4ExcitedStringDecay(new G4QGSMFragmentation);
   theStringModel->SetFragmentationModel(theStringDecay);

   theCascade = new G4GeneratorPrecompoundInterface;
   thePreEquilib = new G4PreCompoundModel(new G4ExcitationHandler);
   theCascade->SetDeExcitation(thePreEquilib);  

   theModel->SetTransport(theCascade);
   theModel->SetHighEnergyGenerator(theStringModel);
 }

 void G4QGSPProtonBuilder::
 Build(G4ProtonInelasticProcess * aP)
 {
// G4cout << "adding inelastic Proton in QGSP" << G4endl;
   aP->AddDataSet(&theXSec);  
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(100*TeV);
   aP->RegisterMe(theModel);
 }

 void G4QGSPProtonBuilder::
 Build(G4HadronElasticProcess * )
 {
 }

 G4QGSPProtonBuilder::
 ~G4QGSPProtonBuilder() 
 {
   delete thePreEquilib;
   delete theCascade;
   delete theStringDecay;
   delete theStringModel;
 }

 // 2002 by J.P. Wellisch
