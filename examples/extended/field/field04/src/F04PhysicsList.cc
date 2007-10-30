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
//

#include "F04PhysicsList.hh"
#include "F04PhysicsListMessenger.hh"

#include "F04ExtraPhysics.hh"
//#include "F04OpticalPhysics.hh"

#include "G4DecayPhysics.hh"

#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics.hh"

#include "G4EmExtraPhysics.hh"
#include "G4EmProcessOptions.hh"

#include "G4HadronElasticPhysics.hh"
#include "G4HadronQElasticPhysics.hh"
#include "G4HadronHElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4QStoppingPhysics.hh"
#include "G4LHEPStoppingPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4IonPhysics.hh"

#include "HadronPhysicsFTFP.hh"
#include "HadronPhysicsFTFC.hh"

#include "HadronPhysicsLHEP.hh"
#include "HadronPhysicsLHEP_BERT.hh"
#include "HadronPhysicsLHEP_EMV.hh"
#include "HadronPhysicsLHEP_PRECO_HP.hh"

#include "HadronPhysicsQGSC.hh"
#include "HadronPhysicsQGSC_EFLOW.hh"

#include "HadronPhysicsQGSP.hh"
#include "HadronPhysicsQGSP_BERT.hh"
#include "HadronPhysicsQGSP_BERT_HP.hh"
#include "HadronPhysicsQGSP_BERT_TRV.hh"
#include "HadronPhysicsQGSP_BIC.hh"
#include "HadronPhysicsQGSP_BIC_HP.hh"

#include "G4HadronInelasticQBBC.hh"
#include "G4HadronInelasticQLHEP.hh"

#include "G4IonPhysics.hh"

#include "G4HadronProcessStore.hh"

#include "G4LossTableManager.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "F04StepMax.hh"

#include "G4ProcessTable.hh"

#include "G4PionDecayMakeSpin.hh"
#include "G4DecayWithSpin.hh"

#include "G4DecayTable.hh"
#include "G4MuonDecayChannelWithSpin.hh"
#include "G4MuonRadiativeDecayChannelWithSpin.hh"

F04PhysicsList::F04PhysicsList(G4String physicsList) : G4VModularPhysicsList() 
{
    G4LossTableManager::Instance();

    defaultCutValue  = 1.*mm;
    fCutForGamma     = defaultCutValue;
    fCutForElectron  = defaultCutValue;
    fCutForPositron  = defaultCutValue;
    fDump            = false;

    fMessenger = new F04PhysicsListMessenger(this);

    fEMPhysics     = new PhysicsListVector();
    fHadronPhysics = new PhysicsListVector();

    // Particles
    fParticleList = new G4DecayPhysics("decays");

    // EM physics
    AddPhysicsList("emstandard_opt1");

    // Set the default hadronic physics.
    AddPhysicsList(physicsList);

    stepMaxProcess = new F04StepMax();
}

F04PhysicsList::~F04PhysicsList()
{
    delete fMessenger;
    delete fParticleList;

    ClearEMPhysics();
    ClearHadronPhysics();

    delete fEMPhysics;
    delete fHadronPhysics;

    delete stepMaxProcess;
}

void F04PhysicsList::ClearEMPhysics()
{
    for (PhysicsListVector::iterator p  = fEMPhysics->begin();
                                     p != fEMPhysics->end(); ++p) {
        delete (*p);
    }
    fEMPhysics->clear();
}

void F04PhysicsList::ClearHadronPhysics()
{
    for (PhysicsListVector::iterator p  = fHadronPhysics->begin();
                                     p != fHadronPhysics->end(); ++p) {
        delete (*p);
    }
    fHadronPhysics->clear();
}

void F04PhysicsList::ConstructParticle()
{
    fParticleList->ConstructParticle();

    G4DecayTable* MuonPlusDecayTable = new G4DecayTable();
    MuonPlusDecayTable -> Insert(new
                           G4MuonDecayChannelWithSpin("mu+",0.986));
    MuonPlusDecayTable -> Insert(new
                           G4MuonRadiativeDecayChannelWithSpin("mu+",0.014));
    G4MuonPlus::MuonPlusDefinition() -> SetDecayTable(MuonPlusDecayTable);

    G4DecayTable* MuonMinusDecayTable = new G4DecayTable();
    MuonMinusDecayTable -> Insert(new
                            G4MuonDecayChannelWithSpin("mu-",0.986));
    MuonMinusDecayTable -> Insert(new
                            G4MuonRadiativeDecayChannelWithSpin("mu-",0.014));
    G4MuonMinus::MuonMinusDefinition() -> SetDecayTable(MuonMinusDecayTable);

}

void F04PhysicsList::ConstructProcess()
{
    AddTransportation();

    for (PhysicsListVector::iterator p  = fEMPhysics->begin();
                                     p != fEMPhysics->end(); ++p) {
        (*p)->ConstructProcess();
    }

    fParticleList->ConstructProcess();

    G4DecayWithSpin* decayWithSpin = new G4DecayWithSpin();

    G4ProcessTable* processTable = G4ProcessTable::GetProcessTable();

    G4VProcess* decay;
    decay = processTable->FindProcess("Decay",G4MuonPlus::MuonPlus());

    G4ProcessManager* fManager;
    fManager = G4MuonPlus::MuonPlus()->GetProcessManager();

    if (fManager) {
      if (decay) fManager->RemoveProcess(decay);
      fManager->AddProcess(decayWithSpin);
      // set ordering for PostStepDoIt and AtRestDoIt
      fManager ->SetProcessOrdering(decayWithSpin, idxPostStep);
      fManager ->SetProcessOrdering(decayWithSpin, idxAtRest);
    }

    decay = processTable->FindProcess("Decay",G4MuonMinus::MuonMinus());

    fManager = G4MuonMinus::MuonMinus()->GetProcessManager();

    if (fManager) {
      if (decay) fManager->RemoveProcess(decay);
      fManager->AddProcess(decayWithSpin);
      // set ordering for PostStepDoIt and AtRestDoIt
      fManager ->SetProcessOrdering(decayWithSpin, idxPostStep);
      fManager ->SetProcessOrdering(decayWithSpin, idxAtRest);
    }

    G4PionDecayMakeSpin* poldecay = new G4PionDecayMakeSpin();

    decay = processTable->FindProcess("Decay",G4PionPlus::PionPlus());

    fManager = G4PionPlus::PionPlus()->GetProcessManager();

    if (fManager) {
      if (decay) fManager->RemoveProcess(decay);
      fManager->AddProcess(poldecay);
      // set ordering for PostStepDoIt and AtRestDoIt
      fManager ->SetProcessOrdering(poldecay, idxPostStep);
      fManager ->SetProcessOrdering(poldecay, idxAtRest);
    }

    decay = processTable->FindProcess("Decay",G4PionMinus::PionMinus());

    fManager = G4PionMinus::PionMinus()->GetProcessManager();

    if (fManager) {
      if (decay) fManager->RemoveProcess(decay);
      fManager->AddProcess(poldecay);
      // set ordering for PostStepDoIt and AtRestDoIt
      fManager ->SetProcessOrdering(poldecay, idxPostStep);
      fManager ->SetProcessOrdering(poldecay, idxAtRest);
    }

    for (PhysicsListVector::iterator p  = fHadronPhysics->begin();
                                     p != fHadronPhysics->end();++p) {
        (*p)->ConstructProcess();
    }

    AddStepMax();

    if (fDump) G4HadronProcessStore::Instance()->Dump(1);
}

void F04PhysicsList::RemoveFromEMPhysicsList(const G4String& name)
{
    G4bool success = false;
    for (PhysicsListVector::iterator p  = fEMPhysics->begin();
                                     p != fEMPhysics->end(); ++p) {
        G4VPhysicsConstructor* e = (*p);
        if (e->GetPhysicsName() == name) {
           fEMPhysics->erase(p);
           success = true;
           break;
        }
    }
    if (!success) {
       std::ostringstream message;
       message << "PhysicsList::RemoveFromEMPhysicsList "<< name << "not found";
       G4Exception(message.str().c_str());
    }
}

void F04PhysicsList::RemoveFromHadronPhysicsList(const G4String& name)
{
    G4bool success = false;
    for (PhysicsListVector::iterator p  = fHadronPhysics->begin();
                                     p != fHadronPhysics->end(); ++p) {
        G4VPhysicsConstructor* e = (*p);
        if (e->GetPhysicsName() == name) {
           fHadronPhysics->erase(p);
           success = true;
           break;
        }
    }
    if (!success) {
       G4String message("F04PhysicsList::RemoveFromHadronPhysicsList \"");
       message += name;
       message += " not found \"";
       G4Exception(message);
    }
}

void F04PhysicsList::AddPhysicsList(const G4String& name)
{
    if (verboseLevel>0) {
        G4cout << "F04PhysicsList::AddPhysicsList: <" << name 
               << ">" << G4endl;
    }

    if (name == "emstandard_opt1") {

       ClearEMPhysics();
       fEMPhysics->push_back(new G4EmStandardPhysics_option1());
       fEMPhysics->push_back(new F04ExtraPhysics());
//       fEMPhysics->push_back(new F04OpticalPhysics());

    } else if (name == "emstandard_opt2") {

       ClearEMPhysics();
       fEMPhysics->push_back(new G4EmStandardPhysics_option2());
       fEMPhysics->push_back(new F04ExtraPhysics());
//       fEMPhysics->push_back(new F04OpticalPhysics());

   } else if (name == "FTFC") {

       SetStandardList(false, false);
       fHadronPhysics->push_back(new HadronPhysicsFTFC("hadron",true));

    } else if (name == "FTFP") {

       SetStandardList(false, false);
       fHadronPhysics->push_back(new HadronPhysicsFTFP("hadron",true));

    } else if (name == "FTFP_EMV") {

       AddPhysicsList("emstandard_opt1");
       AddPhysicsList("FTFP");

    } else if (name == "LHEP") {

       ClearHadronPhysics();
       fHadronPhysics->push_back(new G4EmExtraPhysics("gamma_nuc"));
       fHadronPhysics->push_back(new G4HadronElasticPhysics("LElastic",
                                                           verboseLevel,
                                                           false));
       fHadronPhysics->push_back(new HadronPhysicsLHEP("hadron"));
       fHadronPhysics->push_back(new G4IonPhysics("ion"));
       fDump = true;
 
    } else if (name == "LHEP_BERT") {

       ClearHadronPhysics();
       fHadronPhysics->push_back(new G4EmExtraPhysics("gamma_nuc"));
       fHadronPhysics->push_back(new G4HadronElasticPhysics("LElastic",
                                                           verboseLevel,
                                                           false));
       fHadronPhysics->push_back(new HadronPhysicsLHEP_BERT("hadron"));
       fHadronPhysics->push_back(new G4QStoppingPhysics("stopping",
                                                       verboseLevel));
       fHadronPhysics->push_back(new G4IonPhysics("ion"));
       fDump = true;

    } else if (name == "LHEP_EMV") {

       ClearHadronPhysics();
       AddPhysicsList("emstandard_opt1");
       fHadronPhysics->push_back(new G4EmExtraPhysics("gamma_nuc"));
       fHadronPhysics->push_back(new G4HadronElasticPhysics("LElastic",
                                                           verboseLevel,
                                                           false));
       fHadronPhysics->push_back(new HadronPhysicsLHEP_EMV("hadron"));
       fHadronPhysics->push_back(new G4QStoppingPhysics("stopping",
							verboseLevel));
       fHadronPhysics->push_back(new G4IonPhysics("ion"));
       fDump = true;

    } else if (name == "LHEP_PRECO_HP") {

       ClearHadronPhysics();
       fHadronPhysics->push_back(new G4EmExtraPhysics("gamma_nuc"));
       fHadronPhysics->push_back(new G4HadronElasticPhysics("LElastic",
                                                           verboseLevel,
                                                           false));
       fHadronPhysics->push_back(new HadronPhysicsLHEP_PRECO_HP("hadron"));
       fHadronPhysics->push_back(new G4QStoppingPhysics("stopping",
                                                       verboseLevel));
       fHadronPhysics->push_back(new G4IonPhysics("ion"));
       fDump = true;

    } else if (name == "QBBC") {

       SetStandardList(false,true);
       fHadronPhysics->push_back(new G4HadronInelasticQBBC("QBBC",
                                                          verboseLevel,
                                                          false,false,
                                                          false,false,true));

    } else if (name == "QBBCG") {

       SetStandardList(false, false);
       fHadronPhysics->push_back(new G4HadronInelasticQBBC("QBBC",
                                                          verboseLevel,
                                                          false,false,
                                                          false,false,false));

    } else if (name == "QBEC") {

       SetStandardList(false,true);
       fHadronPhysics->push_back(new G4HadronInelasticQBBC("QBBC",
                                                          verboseLevel,
                                                          false,true,
                                                          false,false,true));

    } else if (name == "QBBC_HP") {

       SetStandardList(true,true);
       fHadronPhysics->push_back(new G4HadronInelasticQBBC("QBBC",
                                                          verboseLevel,
                                                          false,false,
                                                          false,true,true));

    } else if (name == "QBEC_HP") {

       SetStandardList(true,true);
       fHadronPhysics->push_back(new G4HadronInelasticQBBC("QBBC",
                                                          verboseLevel,
                                                          false,true,
                                                          false,true,true));

    } else if (name == "QGSC") {

       ClearHadronPhysics();
       fHadronPhysics->push_back(new G4EmExtraPhysics("gamma_nuc"));
       fHadronPhysics->push_back(new G4HadronQElasticPhysics("QElastic",
                                                            verboseLevel));
       fHadronPhysics->push_back(new HadronPhysicsQGSC("hadron",true));
       fHadronPhysics->push_back(new G4QStoppingPhysics("stopping",
                                                       verboseLevel,false));
       fHadronPhysics->push_back(new G4IonPhysics("ion"));
       fHadronPhysics->push_back(new G4NeutronTrackingCut());
       fDump = true;

    } else if (name == "QGSC_EFLOW") {

       ClearHadronPhysics();
       fHadronPhysics->push_back(new G4EmExtraPhysics("gamma_nuc"));
       fHadronPhysics->push_back(new G4HadronQElasticPhysics("QElastic",
                                                            verboseLevel));
       fHadronPhysics->push_back(new HadronPhysicsQGSC_EFLOW("hadron",true));
       fHadronPhysics->push_back(new G4QStoppingPhysics("stopping",
                                                       verboseLevel,false));
       fHadronPhysics->push_back(new G4IonPhysics("ion"));
       fHadronPhysics->push_back(new G4NeutronTrackingCut());
       fDump = true;

    } else if (name == "QGSC_EMV") {

       AddPhysicsList("emstandard_opt1");
       AddPhysicsList("QGSC");

    } else if (name == "QGSP") {

       SetStandardList(false, false);
       fHadronPhysics->push_back(new HadronPhysicsQGSP("hadron",true));

    } else if (name == "QGSP_BERT") {

       SetStandardList(false, false);
       fHadronPhysics->push_back(new HadronPhysicsQGSP_BERT("hadron",true));

    } else if (name == "QGSP_BERT_EMV") {

       AddPhysicsList("emstandard_opt1");
       AddPhysicsList("QGSP_BERT");

    } else if (name == "QGSP_BERT_HP") {

       SetStandardList(true, false);
       fHadronPhysics->push_back(new HadronPhysicsQGSP_BERT_HP("hadron",true));

    } else if (name == "QGSP_BERT_NQE") {

       SetStandardList(false, false);
       fHadronPhysics->push_back(new HadronPhysicsQGSP_BERT("hadron",false));

    } else if (name == "QGSP_BERT_TRV") {

       SetStandardList(false, false);
       fHadronPhysics->push_back(new HadronPhysicsQGSP_BERT_TRV("hadron",true));

    } else if (name == "QGSP_BIC") {

       SetStandardList(false, false);
       fHadronPhysics->push_back(new HadronPhysicsQGSP_BIC("hadron",true));

    } else if (name == "QGSP_BIC_HP") {

       SetStandardList(true, false);
       fHadronPhysics->push_back(new HadronPhysicsQGSP_BIC_HP("hadron",true));

    } else if (name == "QGSP_EMV") {

       AddPhysicsList("emstandard_opt1");
       AddPhysicsList("QGSP");

    } else if (name == "QGSP_EMV_NQE") {

       AddPhysicsList("emstandard_opt1");
       AddPhysicsList("QGSP_NQE");

    } else if (name == "QGSP_EMX") {

       AddPhysicsList("emstandard_opt2");
       AddPhysicsList("QGSP");

    } else if (name == "QGSP_NQE") {

       SetStandardList(false, false);
       fHadronPhysics->push_back(new HadronPhysicsQGSP("hadron",false));

    } else if (name == "QGSP_QEL") {

       ClearHadronPhysics();
       fHadronPhysics->push_back(new G4EmExtraPhysics("gamma_nuc"));
       fHadronPhysics->push_back(new G4HadronQElasticPhysics("QElastic",
                                                            verboseLevel));
       fHadronPhysics->push_back(new HadronPhysicsQGSP("hadron",true));
       fHadronPhysics->push_back(new G4QStoppingPhysics("stopping",
                                                       verboseLevel));
       fHadronPhysics->push_back(new G4IonPhysics("ion"));
       fHadronPhysics->push_back(new G4NeutronTrackingCut());
       fDump = true;

    } else {

       G4cout << "F04PhysicsList::AddPhysicsList: <" << name << ">"
              << " is not defined" << G4endl;
    }
}

void F04PhysicsList::SetStandardList(G4bool flagHP, G4bool glauber)
{
    ClearHadronPhysics();

    fHadronPhysics->push_back(new G4EmExtraPhysics("gamma_nuc"));
    fHadronPhysics->push_back(new G4HadronElasticPhysics("elastic",
                                                        verboseLevel,
                                                        flagHP,glauber));
    fHadronPhysics->push_back(new  G4QStoppingPhysics("stopping",
                                                     verboseLevel));
    fHadronPhysics->push_back(new G4IonBinaryCascadePhysics("binary_ion"));
    fHadronPhysics->push_back(new G4NeutronTrackingCut());

    fDump = true;
}

void F04PhysicsList::SetCuts()
{
    if (verboseLevel >0) {
        G4cout << "F04PhysicsList::SetCuts:";
        G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length")
               << G4endl;
    }

    // set cut values for gamma at first and for e- second and next for e+,
    // because some processes for e+/e- need cut values for gamma
    SetCutValue(fCutForGamma, "gamma");
    SetCutValue(fCutForElectron, "e-");
    SetCutValue(fCutForPositron, "e+");

    if (verboseLevel>0) DumpCutValuesTable();
}

void F04PhysicsList::SetCutForGamma(G4double cut)
{
    fCutForGamma = cut;
    SetParticleCuts(fCutForGamma, G4Gamma::Gamma());
}

void F04PhysicsList::SetCutForElectron(G4double cut)
{
    fCutForElectron = cut;
    SetParticleCuts(fCutForElectron, G4Electron::Electron());
}

void F04PhysicsList::SetCutForPositron(G4double cut)
{
    fCutForPositron = cut;
    SetParticleCuts(fCutForPositron, G4Positron::Positron());
}

void F04PhysicsList::List()
{
    G4cout << "### PhysicsLists available: FTFC FTFP FTFP_EMV"
           << " LHEP LHEP_BERT LHEP_EMV"
           << G4endl;
    G4cout << "                            LHEP_PRECO_HP QBBC QBBCG"
           << " QBEC QBBC_HP QBEC_HP"
           << G4endl;
    G4cout << "                            QGSC QGSC_EFLOW QGSC_EMV"
           << " QGSP QGSP_BERT QGSP_BER_EMV"
           << G4endl;
    G4cout << "                            QGSP_BERT_HP QGSP_BERT_NQE"
           << " QGSP_BERT_TRV"
           << G4endl;
    G4cout << "                            QGSP_BIC QGSP_BIC_HP QGSP_EMV"
           << " QGSP_EMV_NQE"
           << G4endl;
    G4cout << "                            QGSC_EMX QGSP_NQE QGSP_QEL"
           << G4endl;
}

void F04PhysicsList::SetStepMax(G4double step)
{
  MaxChargedStep = step ;
  stepMaxProcess->SetStepMax(MaxChargedStep);
}

F04StepMax* F04PhysicsList::GetStepMaxProcess()
{
  return stepMaxProcess;
}

void F04PhysicsList::AddStepMax()
{
  // Step limitation seen as a process

  theParticleIterator->reset();
  while ((*theParticleIterator)()){
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();

      if (stepMaxProcess->IsApplicable(*particle) && !particle->IsShortLived())
      {
         if (pmanager) pmanager ->AddDiscreteProcess(stepMaxProcess);
      }
  }
}
