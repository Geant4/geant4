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
// $Id: HadrontherapyProtonElastic.cc; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// (b) INFN Section of Genova, Genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "HadrontherapyProtonElastic.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LElastic.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4PiMinusAbsorptionAtRest.hh"
#include "G4PiNuclearCrossSection.hh" 
#include "G4LFission.hh"
#include "G4LCapture.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4LETritonInelastic.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4LEAlphaInelastic.hh"
#include "G4HadronElasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4BinaryLightIonReaction.hh" 
#include "G4HadronInelasticProcess.hh"
//
// ELASTIC SCATTERING ACTIVE ONLY
//
HadrontherapyProtonElastic::HadrontherapyProtonElastic(const G4String& name): 
  G4VPhysicsConstructor(name)
{
  G4cout<<"****** Proton Elastic Physics List is active !!!!!! ******"
	<<G4endl;

}

HadrontherapyProtonElastic::~HadrontherapyProtonElastic()
{}

void HadrontherapyProtonElastic::ConstructProcess()
{
 G4ParticleDefinition* particle = 0;
 G4ProcessManager* pmanager = 0;

// ELASTIC SCATTERING
 
 //  ELASTIC SCATTERING 
  // FOR PROTON, NEUTRON, IONS 
  G4LElastic* elasticScattering_model = new G4LElastic();
  G4HadronElasticProcess* elastic_scattering = new G4HadronElasticProcess();
  elastic_scattering -> RegisterMe(elasticScattering_model);

// ELASTIC SCATTERING

  ////////////
  // Proton //
  ////////////

  particle = G4Proton::Proton();
  pmanager = particle -> GetProcessManager();
  
  // Activate the proton elastic scattering 
  pmanager -> AddDiscreteProcess(elastic_scattering); 

  /////////////
  // Neutron //
  /////////////

  particle = G4Neutron::Neutron();
  pmanager = particle -> GetProcessManager();
  pmanager -> AddDiscreteProcess(elastic_scattering); 
  
  ////////////////////
  // HADRON CAPTURE //
  ////////////////////

//   G4HadronCaptureProcess* neutronCapture = new G4HadronCaptureProcess();

//   G4LCapture* capture_model = new G4LCapture();

//   // Set the energy range for the capture model
//   capture_model -> SetMinEnergy(neutronLowLimit);
//   capture_model -> SetMaxEnergy(neutronHighLimit);

//   // Register the capture model
//   neutronCapture -> RegisterMe(capture_model);
//   // Active the capture process
//   pmanager -> AddDiscreteProcess(neutronCapture);

 ////////////// 
 // FISSION  //
 /////////////

//   // Process for induced fission
//   G4HadronFissionProcess* fission = new G4HadronFissionProcess();
//   //Final state production model for induced fission
//   G4LFission* fission_model = new G4LFission();
//   // Set the energy range for the fission model
//   fission_model -> SetMinEnergy(neutronLowLimit);
//   fission_model -> SetMaxEnergy(neutronHighLimit);
//   // Register the fission model
//   fission -> RegisterMe(fission_model); 
//   // Active the fission process
//   pmanager -> AddDiscreteProcess(fission);  

  ////////////////		  
  // Pions plus //
  ////////////////
  particle = G4PionPlus::PionPlus(); 
  pmanager = particle -> GetProcessManager();
  pmanager -> AddDiscreteProcess(elastic_scattering);

  ////////////////
  // Pion Minus //
  ///////////////

  particle = G4PionMinus::PionMinus();
  pmanager = particle -> GetProcessManager();
  pmanager -> AddDiscreteProcess(elastic_scattering); 
   
  ///////////////
  // Deuteron //
  //////////////

  particle = G4Deuteron::Deuteron();
  pmanager = particle -> GetProcessManager();

  // Active the Hadron Elastic Process
  pmanager -> AddDiscreteProcess(elastic_scattering); 

  ////////////
  // Triton //
  ////////////
  particle = G4Triton::Triton();
  pmanager = particle -> GetProcessManager();
  // Active the triton elastic scattering process
  pmanager -> AddDiscreteProcess(elastic_scattering);

  ///////////
  // Alpha //
  //////////
  particle = G4Alpha::Alpha();
  pmanager = particle -> GetProcessManager();
  // Active the alpha elastic scattering
  pmanager -> AddDiscreteProcess(elastic_scattering); 

  // He3
  // particle = G4He3::He3();

  // G4HadronInelasticProcess* He3inelasticProcess = 
  //   new G4HadronInelasticProcess("He3Inelastic",particle);
  
  //   G4BinaryLightIonReaction * ionBinaryCascade= new G4BinaryLightIonReaction;
 
  // He3inelasticProcess -> AddDataSet(TripathiCrossSection);
  //He3inelasticProcess -> AddDataSet(aShen);
  //He3inelasticProcess -> RegisterMe(ionBinaryCascade);
 
  //pmanager = particle -> GetProcessManager();
  //pmanager -> AddDiscreteProcess(He3inelasticProcess);
  //pmanager -> AddDiscreteProcess(elastic_scattering); 

}



