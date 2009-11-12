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
//    **************************************
//    *                                    *
//    *    RemSimHadronicBinary.cc        *
//    *                                    *
//    **************************************
//
// $Id: RemSimHadronicBinary.cc,v 1.8 2009-11-12 05:12:18 cirrone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author : Susanna Guatelli, guatelli@ge.infn.it

// Code review: MGP, 7 November 2006 (to be completed)
// 
#include "RemSimHadronicBinary.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LElastic.hh"
//#include "G4CascadeInterface.hh"
//#include "G4PreCompoundModel.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEAlphaInelastic.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"
#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4BinaryCascade.hh"

//
// BINARY PHYSICS LIST
//

RemSimHadronicBinary::RemSimHadronicBinary(const G4String& name): 
G4VPhysicsConstructor(name)
{}

RemSimHadronicBinary::~RemSimHadronicBinary()
{}

void RemSimHadronicBinary::ConstructProcess()
{
  // Physics for proton, neutron, pion+ and pion-
  
  // Elastic scattering: LElastic model
  G4LElastic* elasticModel = new G4LElastic();
  G4HadronElasticProcess* elasticScattering = new G4HadronElasticProcess();
  elasticScattering->RegisterMe(elasticModel);

  // Inelastic scattering: Binary model up to  10. GeV 
  G4BinaryCascade* binaryModel = new G4BinaryCascade();
  // Energy limit of the Binary model
  G4double binaryHighEnergyLimit = 10. * GeV;
  binaryModel->SetMaxEnergy(binaryHighEnergyLimit);

  // Inelastic scattering: LEP model between 8. * GeV and 25. * GeV 
  G4LEProtonInelastic* LEPProtonModel = new G4LEProtonInelastic();
  G4LENeutronInelastic* LEPNeutronModel = new G4LENeutronInelastic();
  G4LEPionPlusInelastic* LEPPionPlusModel = new G4LEPionPlusInelastic();
  G4LEPionMinusInelastic* LEPPionMinusModel = new G4LEPionMinusInelastic();
  // Set the energy limits
  G4double LEPLowEnergyLimit = 8. * GeV;
  G4double LEPHighEnergyLimit = 25. * GeV;
  LEPProtonModel->SetMinEnergy(LEPLowEnergyLimit);
  LEPProtonModel->SetMaxEnergy(LEPHighEnergyLimit);
  LEPNeutronModel->SetMinEnergy(LEPLowEnergyLimit);
  LEPNeutronModel->SetMaxEnergy(LEPHighEnergyLimit);

  //  no intranuclear transport activated for pions
  // at this stage; further tests on Binary Cascade for pions needed 
 
  G4double LEPPionLowEnergyLimit = 0. * MeV;
  LEPPionPlusModel->SetMinEnergy(LEPPionLowEnergyLimit);
  LEPPionPlusModel->SetMaxEnergy(LEPHighEnergyLimit);
  LEPPionMinusModel->SetMinEnergy(LEPPionLowEnergyLimit);
  LEPPionMinusModel->SetMaxEnergy(LEPHighEnergyLimit);

  // Inelastic scattering: QGSP model between 20 GeV and 100 TeV
  QGSPModel = new G4TheoFSGenerator();
  // Set the energy limits
  G4double QGSPLowEnergyLimit = 20.* GeV;
  G4double QGSPHighEnergyLimit = 100.* GeV;
  QGSPModel->SetMinEnergy(QGSPLowEnergyLimit);
  QGSPModel->SetMaxEnergy(QGSPHighEnergyLimit);

  theCascade = new G4GeneratorPrecompoundInterface();
  thePreEquilib = new G4PreCompoundModel(&theHandler);
  theCascade->SetDeExcitation(thePreEquilib);
  QGSPModel->SetTransport(theCascade);

  // Set the fragmentation 
  theFragmentation = new G4QGSMFragmentation();
  theStringDecay = new G4ExcitedStringDecay(theFragmentation);
  theStringModel.SetFragmentationModel(theStringDecay);
  QGSPModel->SetHighEnergyGenerator(&theStringModel);

  // ---------------------------------------------------------------------------------------------
  // Proton processes
  G4ParticleDefinition* proton = G4Proton::ProtonDefinition();
  G4ProcessManager* protonProcessManager = proton->GetProcessManager();
 
  // Proton inelastic scattering
  G4ProtonInelasticProcess* protonInelasticProcess = new G4ProtonInelasticProcess();
  // Activate the cross-sections for proton nuclear scattering up to 20 GeV
  //G4ProtonInelasticCrossSection protonCrossSection;
  protonInelasticProcess->AddDataSet(&protonCrossSection);
  // Set the models
  protonInelasticProcess->RegisterMe(binaryModel);
  protonInelasticProcess->RegisterMe(LEPProtonModel);
  
  protonInelasticProcess->RegisterMe(QGSPModel);
  // Activate the inelastic scattering
  protonProcessManager->AddDiscreteProcess(protonInelasticProcess);
  // Activate the elastic scattering
  protonProcessManager->AddDiscreteProcess(elasticScattering);

  //------------------------------------------------------
  // Pion Plus processes
  G4ParticleDefinition* piPlus = G4PionPlus::PionPlusDefinition();
  G4ProcessManager* pionPlusProcessManager = piPlus->GetProcessManager();
 
  // Define the inelastic scattering for pion plus 
  G4PionPlusInelasticProcess* pionPlusInelasticProcess = new G4PionPlusInelasticProcess();
  // Set the cross section
 // G4PiNuclearCrossSection pionCrossSection;
  pionPlusInelasticProcess->AddDataSet(&pionCrossSection);
  // Register the models
  pionPlusInelasticProcess->RegisterMe(LEPPionPlusModel);
  pionPlusInelasticProcess->RegisterMe(QGSPModel);
  // Activate the inelastic scattering
  pionPlusProcessManager->AddDiscreteProcess(pionPlusInelasticProcess);
  // Activate the elastic process
  pionPlusProcessManager->AddDiscreteProcess(elasticScattering);

  //------------------------------------------------------------
  // Pion Minus processes
  G4ParticleDefinition* piMinus = G4PionMinus::PionMinusDefinition();
  G4ProcessManager* pionMinusProcessManager = piMinus->GetProcessManager();
  
  // Define the inelastic processes for pion minus
  G4PionMinusInelasticProcess* pionMinusInelasticProcess = new G4PionMinusInelasticProcess();
  // Set the cross section
  pionMinusInelasticProcess->AddDataSet(&pionCrossSection);
  // Register the models
  pionMinusInelasticProcess->RegisterMe(LEPPionMinusModel);
  pionMinusInelasticProcess->RegisterMe(QGSPModel);
  // Activate the inelastic scattering
  pionMinusProcessManager->AddDiscreteProcess(pionMinusInelasticProcess);
  // Activate the elastic scattering
  pionMinusProcessManager->AddDiscreteProcess(elasticScattering);
 
  //-----------------------------------------------------
  // Neutron processes
  G4ParticleDefinition* neutron = G4Neutron::NeutronDefinition();
  G4ProcessManager* neutronProcessManager = neutron->GetProcessManager();

  // Inelastic process
  G4NeutronInelasticProcess* neutronInelasticProcess = new G4NeutronInelasticProcess();
  // Set the cross section
  //G4NeutronInelasticCrossSection neutronCrossSection;
  neutronInelasticProcess->AddDataSet(&neutronCrossSection);
  // Set the models
  neutronInelasticProcess->RegisterMe(binaryModel);
  neutronInelasticProcess->RegisterMe(LEPNeutronModel);
  neutronInelasticProcess->RegisterMe(QGSPModel);
 // Activate the neutron inelastic scattering
  neutronProcessManager->AddDiscreteProcess(neutronInelasticProcess);
  // Activate the neutron elastic scattering
  neutronProcessManager->AddDiscreteProcess(elasticScattering);

  // Neutron capture process
  G4HadronCaptureProcess* neutronCapture = new G4HadronCaptureProcess();
  // Final state production model for capture of neutral hadrons in nuclei
  G4LCapture* captureModel = new G4LCapture();
  // Set the energy range for the capture model
  G4double neutronLowEnergyLimit = 0. * MeV;
  G4double neutronHighEnergyLimit = 100. * TeV;
  captureModel->SetMinEnergy(neutronLowEnergyLimit);
  captureModel->SetMaxEnergy(neutronHighEnergyLimit);
  // Activate the neutron capture model
  neutronCapture->RegisterMe(captureModel);
  // Activate the neutron capture process
  neutronProcessManager->AddDiscreteProcess(neutronCapture);

  // Process for induced fission
   G4HadronFissionProcess* fission = new G4HadronFissionProcess();
  //Final state production model for induced fission
  G4LFission* fissionModel = new G4LFission();
  // Set the energy range for the fission model
  fissionModel->SetMinEnergy(neutronLowEnergyLimit);
  fissionModel->SetMaxEnergy(neutronHighEnergyLimit);
  // Register the fission model
  fission->RegisterMe(fissionModel); 
  // Activate the fission process
  neutronProcessManager->AddDiscreteProcess(fission);  

  //--------------------------------------------------------
  // Physics for alpha particles
 
  G4ParticleDefinition* alpha = G4Alpha::AlphaDefinition();
  G4ProcessManager* alphaProcessManager = alpha->GetProcessManager();

  // Cross section data sets
  
  // TRIPATHI CROSS SECTION
  // Implementation of formulas taken from NASA technical paper 3621 by 
  // Tripathi, et al. Cross-sections for ion ion scattering
  G4TripathiCrossSection* tripathi = new G4TripathiCrossSection();

  // IONS SHEN CROSS SECTION
  // Implementation of formulas 
  // Shen et al. Nuc. Phys. A 491 130 (1989) 
  // Total Reaction Cross Section for Heavy-Ion Collisions
  G4IonsShenCrossSection* shen = new G4IonsShenCrossSection();

  G4LEAlphaInelastic* LEPAlphaModel = new G4LEAlphaInelastic();
  // Energy limit of the LEP model for alpha particles
  G4double LEPAlphaHighLimit = 100 * MeV;
  LEPAlphaModel->SetMaxEnergy(LEPAlphaHighLimit);

  G4BinaryLightIonReaction* binaryIonModel = new G4BinaryLightIonReaction();
  // Energy limit of the binary ion model
  G4double binaryIonLowLimit = 80. * MeV;
  G4double binaryIonHighLimit = 400. * GeV; 
  binaryIonModel->SetMinEnergy(binaryIonLowLimit);
  binaryIonModel->SetMaxEnergy(binaryIonHighLimit);

  // Define the alpha inelastic scattering
  G4AlphaInelasticProcess* alphaInelasticProcess = new G4AlphaInelasticProcess();
  // Activate the Tripathi and Shen Cross Section
  alphaInelasticProcess->AddDataSet(tripathi);
  alphaInelasticProcess->AddDataSet(shen);
  // Set the models
  alphaInelasticProcess->RegisterMe(LEPAlphaModel);
  alphaInelasticProcess->RegisterMe(binaryIonModel);
  // Activate the inelastic scattering
  alphaProcessManager->AddDiscreteProcess(alphaInelasticProcess);
  // Activate the elastic scattering
  alphaProcessManager->AddDiscreteProcess(elasticScattering);
}



