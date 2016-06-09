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
//    *    RemSimHadronicPhysics.cc        *
//    *                                    *
//    **************************************
//
// $Id: RemSimHadronicBinary.cc,v 1.6 2006/06/29 16:23:49 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// Author : Susanna Guatelli, guatelli@ge.infn.it
// 
#include "RemSimHadronicBinary.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4ios.hh"

#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4LElastic.hh"
#include "G4CascadeInterface.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEAlphaInelastic.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"
#include "G4TheoFSGenerator.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4BinaryCascade.hh"
#include "G4PiMinusAbsorptionAtRest.hh"

RemSimHadronicBinary::RemSimHadronicBinary(const G4String& name): 
G4VPhysicsConstructor(name)
{
  // Hadronic cross sections

  proton_XC = new G4ProtonInelasticCrossSection();
  neutron_XC = new G4NeutronInelasticCrossSection();
  pion_XC = new G4PiNuclearCrossSection();
 
  // Hadronic cross sections for alpha particles 
  tripathi = new G4TripathiCrossSection();
  shen = new G4IonsShenCrossSection();

  // Hadronic models

  elastic_model = new G4LElastic();

  binarycascade_model = new G4BinaryCascade();
  //bertini_model = new G4CascadeInterface();

  binary_ion_model = new G4BinaryLightIonReaction();
 
  LEP_proton_model = new G4LEProtonInelastic();
  LEP_neutron_model = new G4LENeutronInelastic();
  LEP_pip_model = new G4LEPionPlusInelastic();
  LEP_pim_model = new G4LEPionMinusInelastic();
  LEP_alpha_model = new G4LEAlphaInelastic();
  nfission_model = new G4LFission();
  ncapture_model = new G4LCapture();

  QGSP_model = new G4TheoFSGenerator();
  theCascade = new G4GeneratorPrecompoundInterface();
 
  thePreEquilib = new G4PreCompoundModel(&theHandler);
  theCascade -> SetDeExcitation(thePreEquilib);
  QGSP_model -> SetTransport(theCascade);
  theFragmentation = new G4QGSMFragmentation();
  theStringDecay = new G4ExcitedStringDecay(theFragmentation);
  theStringModel.SetFragmentationModel(theStringDecay);
  QGSP_model -> SetHighEnergyGenerator(&theStringModel);
 
}
RemSimHadronicBinary::~RemSimHadronicBinary()
{}

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"

void RemSimHadronicBinary::ConstructProcess()
{
  // 0- 10. GeV: binary cascade for p, n

  binarycascade_model -> SetMaxEnergy(10.*GeV);

  // 80 MeV - 20 GeV for light ions

  binary_ion_model -> SetMinEnergy(80*MeV);
  binary_ion_model -> SetMaxEnergy(110.*GeV);
 
  LEP_proton_model -> SetMinEnergy(8.*GeV);
  LEP_proton_model -> SetMaxEnergy(25*GeV);
 
  LEP_neutron_model -> SetMinEnergy(8.*GeV);
  LEP_neutron_model -> SetMaxEnergy(25*GeV);
 
  LEP_pip_model -> SetMinEnergy(0.*GeV);
  LEP_pip_model -> SetMaxEnergy(25*GeV);
  LEP_pim_model -> SetMinEnergy(0.*GeV);
  LEP_pim_model -> SetMaxEnergy(25*GeV);
 
  nfission_model -> SetMinEnergy(0*TeV);
  nfission_model -> SetMaxEnergy(100*TeV);
  ncapture_model -> SetMinEnergy(0*TeV);
  ncapture_model -> SetMaxEnergy(100*TeV);

  // Up to 100 MeV for alphas
  LEP_alpha_model -> SetMaxEnergy(100.*MeV);

  // 20 GeV - 100 TeV: QGSP model for p, n, pi+, pi-

  QGSP_model -> SetMinEnergy(20*GeV);
  QGSP_model -> SetMaxEnergy(100*TeV);

  // ******************************************
  // * register models, processes to proton   *
  // ******************************************

  G4ParticleDefinition* proton = G4Proton::ProtonDefinition();
  G4ProcessManager* protMan = proton -> GetProcessManager();

  // Elastic process

  G4HadronElasticProcess* protelProc = new G4HadronElasticProcess();
  protelProc -> RegisterMe(elastic_model);
  protMan -> AddDiscreteProcess(protelProc);

  // Inelastic process

  G4ProtonInelasticProcess* protinelProc = new G4ProtonInelasticProcess();
  protinelProc -> AddDataSet(proton_XC);
  protinelProc -> RegisterMe(binarycascade_model);
  protinelProc -> RegisterMe(LEP_proton_model);
  protinelProc -> RegisterMe(QGSP_model);
  protMan -> AddDiscreteProcess(protinelProc);

  // ******************************************
  // * register models, processes to neutron  *
  // ******************************************

  G4ParticleDefinition* neutron = G4Neutron::NeutronDefinition();
  G4ProcessManager* neutMan = neutron->GetProcessManager();

  // Elastic process

  G4HadronElasticProcess* neutelProc = new G4HadronElasticProcess();
  neutelProc -> RegisterMe(elastic_model);
  neutMan -> AddDiscreteProcess(neutelProc);

  // Inelastic process

  G4NeutronInelasticProcess* neutinelProc = new G4NeutronInelasticProcess();
  neutinelProc -> AddDataSet(neutron_XC);
  neutinelProc -> RegisterMe(binarycascade_model);
  neutinelProc -> RegisterMe(LEP_neutron_model);
  neutinelProc -> RegisterMe(QGSP_model);
  neutMan -> AddDiscreteProcess(neutinelProc);

  G4HadronCaptureProcess* neutcapProc = new G4HadronCaptureProcess();
  neutcapProc -> RegisterMe(ncapture_model);
  neutMan -> AddDiscreteProcess(neutcapProc);

  G4HadronFissionProcess* neutfisProc = new G4HadronFissionProcess();
  neutfisProc -> RegisterMe(nfission_model);
  neutMan -> AddDiscreteProcess(neutfisProc);

  // ******************************************
  // * register models, processes to pi+      *
  // ******************************************

  G4ParticleDefinition* piplus = G4PionPlus::PionPlusDefinition();
  G4ProcessManager* pipMan = piplus -> GetProcessManager();

  // Elastic process

  G4HadronElasticProcess* pipelProc = new G4HadronElasticProcess();
  pipelProc -> RegisterMe(elastic_model);
  pipMan -> AddDiscreteProcess(pipelProc);

  // Inelastic process
    
  G4PionPlusInelasticProcess* pipinelProc = new G4PionPlusInelasticProcess();
  pipinelProc -> AddDataSet(pion_XC);
  pipinelProc -> RegisterMe(LEP_pip_model);
  pipinelProc -> RegisterMe(QGSP_model);                                                
  pipMan -> AddDiscreteProcess(pipinelProc);

  // ******************************************
  // * register models, processes to pi-      *
  // ******************************************

  G4ParticleDefinition* piminus = G4PionMinus::PionMinusDefinition();
  G4ProcessManager* pimMan = piminus -> GetProcessManager();

  // Elastic process

  G4HadronElasticProcess* pimelProc = new G4HadronElasticProcess();
  pimelProc -> RegisterMe(elastic_model);
  pimMan -> AddDiscreteProcess(pimelProc);

  // Inelastic process

  G4PionMinusInelasticProcess* piminelProc = new G4PionMinusInelasticProcess();
  piminelProc -> AddDataSet(pion_XC);
  piminelProc -> RegisterMe(LEP_pim_model);
  piminelProc -> RegisterMe(QGSP_model); 
  //pimMan -> AddRestProcess(new G4PiMinusAbsorptionAtRest, ordDefault);
  pimMan -> AddDiscreteProcess(piminelProc);

  // ******************************************
  // * register models, processes to alpha    *
  // ******************************************

  G4ParticleDefinition* alpha = G4Alpha::AlphaDefinition();
  G4ProcessManager* alfMan = alpha -> GetProcessManager();

  // Elastic process

  G4HadronElasticProcess* alfelProc = new G4HadronElasticProcess();
  alfelProc -> RegisterMe(elastic_model);
  alfMan -> AddDiscreteProcess(alfelProc);

  // Inelastic process

  G4AlphaInelasticProcess* alfinelProc = new G4AlphaInelasticProcess();
  alfinelProc -> AddDataSet(tripathi);
  alfinelProc -> AddDataSet(shen);
  alfinelProc -> RegisterMe(LEP_alpha_model);
  alfinelProc -> RegisterMe(binary_ion_model);
  alfMan -> AddDiscreteProcess(alfinelProc);
}



