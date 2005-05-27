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
//    **************************************
//    *                                    *
//    *    HadrontherapyProtonHadro.cc        *
//    *                                    *
//    **************************************
//
// $Id: HadrontherapyProtonHadro.cc,v 1.5 2005-05-27 15:05:53 cirrone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author : Susanna Guatelli, guatelli@ge.infn.it
// 
#include "HadrontherapyProtonHadro.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4LElastic.hh"
#include "G4BinaryCascade.hh"
#include "G4CascadeInterface.hh"
HadrontherapyProtonHadro::HadrontherapyProtonHadro(const G4String& name): 
G4VPhysicsConstructor(name)
{
}
HadrontherapyProtonHadro::~HadrontherapyProtonHadro()
{}

void HadrontherapyProtonHadro::ConstructProcess()
{
  // ELASTIC SCATTERING
  // ****PROTON, NEUTRON, IONS
  G4LElastic* elastic_Model = new G4LElastic();
  G4HadronElasticProcess* elastic = new G4HadronElasticProcess();
  elastic -> RegisterMe(elastic_Model);

  // INELASTIC SCATTERING
      // Binary Cascade
  G4ParticleDefinition* particle = 0;
  G4ProcessManager* pmanager = 0;
  G4BinaryLightIonReaction* theBC = 0;
  theBC = new G4BinaryLightIonReaction();
  theBC->SetMinEnergy(80*MeV);
  theBC->SetMaxEnergy(20*GeV);

  G4TripathiCrossSection * TripathiCrossSection= new G4TripathiCrossSection;
  G4IonsShenCrossSection * aShen = new G4IonsShenCrossSection;

  // deuteron
  particle = G4Deuteron::Deuteron();
  pmanager = particle->GetProcessManager();
  G4LEDeuteronInelastic* theDIModel = new G4LEDeuteronInelastic;
  theDIModel->SetMaxEnergy(100*MeV);
  theIPdeuteron.AddDataSet(TripathiCrossSection);
  theIPdeuteron.AddDataSet(aShen);
  theIPdeuteron.RegisterMe(theDIModel);
  theIPdeuteron.RegisterMe(theBC);
  pmanager->AddDiscreteProcess(&theIPdeuteron);
  pmanager -> AddDiscreteProcess(elastic); // ELASTIC PROCESSES

  // triton
  particle = G4Triton::Triton();
  pmanager = particle->GetProcessManager();
  G4LETritonInelastic* theTIModel = new G4LETritonInelastic;
  theTIModel->SetMaxEnergy(100*MeV);
  theIPtriton.AddDataSet(TripathiCrossSection);
  theIPtriton.AddDataSet(aShen);
  theIPtriton.RegisterMe(theTIModel);
  theIPtriton.RegisterMe(theBC);
  pmanager->AddDiscreteProcess(&theIPtriton);
  pmanager -> AddDiscreteProcess(elastic); //ELASTIC SCATTERING

  // alpha
  particle = G4Alpha::Alpha();
  pmanager = particle->GetProcessManager();
  G4LEAlphaInelastic* theAIModel = new G4LEAlphaInelastic;
  theAIModel->SetMaxEnergy(100*MeV);
  theIPalpha.AddDataSet(TripathiCrossSection);
  theIPalpha.AddDataSet(aShen);
  theIPalpha.RegisterMe(theAIModel);
  theIPalpha.RegisterMe(theBC);
  pmanager->AddDiscreteProcess(&theIPalpha);
  pmanager -> AddDiscreteProcess(elastic); // ELASTIC SCATTERING

  // Proton PRECOMPOUND

  particle = G4Proton::Proton();
  pmanager = particle->GetProcessManager();
  G4PreCompoundModel* thePreEquilib = new G4PreCompoundModel(&theHandler);
  thePreEquilib->SetMinEnergy(0*MeV);
  thePreEquilib->SetMaxEnergy(10*GeV);
  theIPProton.RegisterMe(thePreEquilib);
  theIPProton.AddDataSet(&thePXSec);
  pmanager->AddDiscreteProcess(&theIPProton);
  pmanager -> AddDiscreteProcess(elastic); // ELASIC SCATTERING

  
  // He3
  particle = G4He3::He3();
  pmanager = particle->GetProcessManager();
  G4BinaryLightIonReaction * theGenIonBC= new G4BinaryLightIonReaction;
  G4HadronInelasticProcess* theIPHe3 =
    new G4HadronInelasticProcess("He3Inelastic",particle);
  theIPHe3->AddDataSet(TripathiCrossSection);
  theIPHe3->AddDataSet(aShen);
  theIPHe3->RegisterMe(theGenIonBC);
  pmanager->AddDiscreteProcess(theIPHe3);
  pmanager -> AddDiscreteProcess(elastic); // ELASTIC SCATTERING
  
  // Neutron
  particle = G4Neutron::Neutron();
  pmanager = particle->GetProcessManager();
  G4PreCompoundModel* thePreEquilib = new G4PreCompoundModel(&theHandler);
  thePreEquilib->SetMinEnergy(0*MeV);
  thePreEquilib->SetMaxEnergy(10*GeV);
  theIPNeutron.RegisterMe(thePreEquilib);
  theIPNeutron.AddDataSet(&theNXSec);
  pmanager -> AddDiscreteProcess(&theIPNeutron);
  pmanager -> AddDiscreteProcess(elastic); // ELASTIC SCATTERING

  //Hadron Capture 
  G4HadronCaptureProcess* neutronCapture = new G4HadronCaptureProcess();
  G4LCapture* capture_model = new G4LCapture();
  capture_model -> SetMinEnergy(0. *TeV);
  capture_model -> SetMaxEnergy(100. *TeV);
  neutronCapture -> RegisterMe(capture_model);
  pmanager -> AddDiscreteProcess(neutronCapture);

  //Fission

  G4HadronFissionProcess* fission = new G4HadronFissionProcess();
  G4LFission* fission_model = new G4LFission();
  fission_model -> SetMinEnergy(0. *TeV);
  fission_model -> SetMaxEnergy(100. *TeV);
  fission -> RegisterMe(fission_model); 
  pmanager -> AddDiscreteProcess(fission);      
}



