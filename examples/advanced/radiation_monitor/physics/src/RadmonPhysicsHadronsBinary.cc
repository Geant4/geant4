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
// File name:     RadmonPhysicsHadronsBinary.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsHadronsBinary.cc,v 1.3 2006/06/29 16:18:47 gunter Exp $
// Tag:           $Name: geant4-09-02 $
//

#include "RadmonPhysicsHadronsBinary.hh"

#include "G4Proton.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4Alpha.hh"

#include "G4ProcessManager.hh"

#include "G4LElastic.hh"
#include "G4TheoFSGenerator.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4QGSMFragmentation.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4HadronElasticProcess.hh"
#include "G4BinaryCascade.hh"
#include "G4LEProtonInelastic.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4LEAlphaInelastic.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"



RadmonVSubPhysicsListWithLabel *                RadmonPhysicsHadronsBinary :: New(void) const
{
 return new RadmonPhysicsHadronsBinary;
}



void                                            RadmonPhysicsHadronsBinary :: ConstructParticle(void)
{
 G4Proton::ProtonDefinition();
 G4PionMinus::PionMinusDefinition();
 G4Alpha::AlphaDefinition();
 G4PionPlus::PionPlusDefinition();
}



void                                            RadmonPhysicsHadronsBinary :: ConstructProcess(void)
{
 // Models
 G4LElastic * elasticModel(new G4LElastic);
 
 G4TheoFSGenerator * qgspModel(new G4TheoFSGenerator());
 G4GeneratorPrecompoundInterface * cascade(new G4GeneratorPrecompoundInterface());
 cascade->SetDeExcitation(new G4PreCompoundModel(new G4ExcitationHandler()));
 qgspModel->SetTransport(cascade);
 G4QGSModel<G4QGSParticipants> * stringModel(new G4QGSModel<G4QGSParticipants>);
 stringModel->SetFragmentationModel(new G4ExcitedStringDecay(new G4QGSMFragmentation()));
 qgspModel->SetHighEnergyGenerator(stringModel);
 qgspModel->SetMinEnergy(20*GeV);
 qgspModel->SetMaxEnergy(100*TeV);
 
 G4PiNuclearCrossSection * piNuclear(new G4PiNuclearCrossSection());


 
 
 
 // Protons
 G4ProcessManager * manager(G4Proton::ProtonDefinition()->GetProcessManager());



 G4HadronElasticProcess * protonElasticProcess(new G4HadronElasticProcess());
 protonElasticProcess->RegisterMe(elasticModel);
 manager->AddDiscreteProcess(protonElasticProcess);



 G4BinaryCascade * protonBinaryCascadeModel(new G4BinaryCascade());
 protonBinaryCascadeModel->SetMaxEnergy(10.*GeV);

 G4LEProtonInelastic * protonLEPModel(new G4LEProtonInelastic());
 protonLEPModel->SetMinEnergy(8.*GeV);
 protonLEPModel->SetMaxEnergy(25.*GeV);

 G4ProtonInelasticProcess * protonInelasticProcess(new G4ProtonInelasticProcess());
 protonInelasticProcess->AddDataSet(new G4ProtonInelasticCrossSection());
 protonInelasticProcess->RegisterMe(protonBinaryCascadeModel);
 protonInelasticProcess->RegisterMe(protonLEPModel);
 protonInelasticProcess->RegisterMe(qgspModel);
 manager->AddDiscreteProcess(protonInelasticProcess);





 // Pi+
 manager=G4PionPlus::PionPlusDefinition()->GetProcessManager();



 G4HadronElasticProcess * pionPlusElasticProcess(new G4HadronElasticProcess());
 pionPlusElasticProcess->RegisterMe(elasticModel);
 manager->AddDiscreteProcess(pionPlusElasticProcess);



 G4LEPionPlusInelastic * pionPlusLEPModel(new G4LEPionPlusInelastic());
 pionPlusLEPModel->SetMinEnergy(0.*GeV);
 pionPlusLEPModel->SetMaxEnergy(25.*GeV);

 G4PionPlusInelasticProcess * pionPlusInelasticProcess(new G4PionPlusInelasticProcess());
 pionPlusInelasticProcess->AddDataSet(piNuclear);
 pionPlusInelasticProcess->RegisterMe(pionPlusLEPModel);
 pionPlusInelasticProcess->RegisterMe(qgspModel);
 manager->AddDiscreteProcess(pionPlusInelasticProcess);





 // Pi-
 manager=G4PionMinus::PionMinusDefinition()->GetProcessManager();



 G4HadronElasticProcess * pionMinusElasticProcess(new G4HadronElasticProcess());
 pionMinusElasticProcess->RegisterMe(elasticModel);
 manager->AddDiscreteProcess(pionMinusElasticProcess);



 G4LEPionMinusInelastic * pionMinusLEPModel(new G4LEPionMinusInelastic());
 pionMinusLEPModel->SetMinEnergy(0.*GeV);
 pionMinusLEPModel->SetMaxEnergy(25.*GeV);

 G4PionMinusInelasticProcess * pionMinusInelasticProcess(new G4PionMinusInelasticProcess());
 pionMinusInelasticProcess->AddDataSet(piNuclear);
 pionMinusInelasticProcess->RegisterMe(pionMinusLEPModel);
 pionMinusInelasticProcess->RegisterMe(qgspModel);
 manager->AddDiscreteProcess(pionMinusInelasticProcess);





 // Alpha
 manager=G4Alpha::AlphaDefinition()->GetProcessManager();



 G4HadronElasticProcess * alphaElasticProcess(new G4HadronElasticProcess());
 alphaElasticProcess->RegisterMe(elasticModel);
 manager->AddDiscreteProcess(alphaElasticProcess);



 G4LEAlphaInelastic * alphaLEPModel(new G4LEAlphaInelastic());
 alphaLEPModel->SetMinEnergy(0.*GeV);
 alphaLEPModel->SetMaxEnergy(25.*GeV);

 G4BinaryLightIonReaction * alphaBinaryModel(new G4BinaryLightIonReaction());
 alphaLEPModel->SetMinEnergy(80.*MeV);
 alphaLEPModel->SetMaxEnergy(110.*GeV);
 
 G4AlphaInelasticProcess * alphaInelasticProcess(new G4AlphaInelasticProcess());
 alphaInelasticProcess->AddDataSet(new G4TripathiCrossSection());
 alphaInelasticProcess->AddDataSet(new G4IonsShenCrossSection());
 alphaInelasticProcess->RegisterMe(alphaLEPModel);
 alphaInelasticProcess->RegisterMe(alphaBinaryModel);
 manager->AddDiscreteProcess(alphaInelasticProcess);
}



void                                            RadmonPhysicsHadronsBinary :: SetCuts(void)
{
}





const RadmonPhysicsInfoList &                   RadmonPhysicsHadronsBinary :: Provides(void) const
{
 if (infoList.GetNPhysicsInfos()==0)
 {
  RadmonPhysicsInfo info;
  
  info.SetProcessName("Elastic");
  info.SetParticleDefinition(G4Proton::ProtonDefinition());
  info.SetMinEnergy(0.*eV);
  info.SetMaxEnergy(DBL_MAX);
  infoList.InsertPhysicsInfo(info);

  info.SetProcessName("Inelastic");
  info.SetParticleDefinition(G4Proton::ProtonDefinition());
  info.SetMinEnergy(0.*eV);
  info.SetMaxEnergy(100.*TeV);
  infoList.InsertPhysicsInfo(info);

  info.SetProcessName("Elastic");
  info.SetParticleDefinition(G4PionPlus::PionPlusDefinition());
  info.SetMinEnergy(0.*eV);
  info.SetMaxEnergy(DBL_MAX);
  infoList.InsertPhysicsInfo(info);

  info.SetProcessName("Inelastic");
  info.SetParticleDefinition(G4PionPlus::PionPlusDefinition());
  info.SetMinEnergy(0.*eV);
  info.SetMaxEnergy(100.*TeV);
  infoList.InsertPhysicsInfo(info);

  info.SetProcessName("Elastic");
  info.SetParticleDefinition(G4PionMinus::PionMinusDefinition());
  info.SetMinEnergy(0.*eV);
  info.SetMaxEnergy(DBL_MAX);
  infoList.InsertPhysicsInfo(info);

  info.SetProcessName("Inelastic");
  info.SetParticleDefinition(G4PionMinus::PionMinusDefinition());
  info.SetMinEnergy(0.*eV);
  info.SetMaxEnergy(100.*TeV);
  infoList.InsertPhysicsInfo(info);

  info.SetProcessName("Elastic");
  info.SetParticleDefinition(G4Alpha::AlphaDefinition());
  info.SetMinEnergy(0.*eV);
  info.SetMaxEnergy(DBL_MAX);
  infoList.InsertPhysicsInfo(info);

  info.SetProcessName("Inelastic");
  info.SetParticleDefinition(G4Alpha::AlphaDefinition());
  info.SetMinEnergy(0.*eV);
  info.SetMaxEnergy(110.*GeV);
  infoList.InsertPhysicsInfo(info);
 }
 
 return infoList;
}
