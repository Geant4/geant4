// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN04HadronPhysics.cc,v 1.1 2001-01-06 06:56:18 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "ExN04HadronPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/iomanip"   


ExN04HadronPhysics::ExN04HadronPhysics(const G4String& name)
                    :  G4VPhysicsConstructor(name)
{
}

ExN04HadronPhysics::~ExN04HadronPhysics()
{
}

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

// Nuclei
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

void ExN04HadronPhysics::ConstructParticle()
{
  //  Construct all mesons
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  //  Construct all barions
  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  

}


#include "G4ProcessManager.hh"


void ExN04HadronPhysics::ConstructProcess()
{
  G4ProcessManager * pManager = 0;
  
  // Elastic Process
  theElasticModel = new G4LElastic();
  theElasticProcess.RegisterMe(theElasticModel);

  // PionPlus
  pManager = G4PionPlus::PionPlus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEPionPlusModel = new G4LEPionPlusInelastic();
  theHEPionPlusModel = new G4HEPionPlusInelastic();
  thePionPlusInelastic.RegisterMe(theLEPionPlusModel);
  thePionPlusInelastic.RegisterMe(theHEPionPlusModel);
  pManager->AddDiscreteProcess(&thePionPlusInelastic);

  pManager->AddProcess(&thePionPlusIonisation, ordInActive,2, 2);

  pManager->AddProcess(&thePionPlusMult);
  pManager->SetProcessOrdering(&thePionPlusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&thePionPlusMult, idxPostStep, 1);

  // PionMinus
  pManager = G4PionMinus::PionMinus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEPionMinusModel = new G4LEPionMinusInelastic();
  theHEPionMinusModel = new G4HEPionMinusInelastic();
  thePionMinusInelastic.RegisterMe(theLEPionMinusModel);
  thePionMinusInelastic.RegisterMe(theHEPionMinusModel);
  pManager->AddDiscreteProcess(&thePionMinusInelastic);

  pManager->AddProcess(&thePionMinusIonisation, ordInActive,2, 2);

  pManager->AddProcess(&thePionMinusMult);
  pManager->SetProcessOrdering(&thePionMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&thePionMinusMult, idxPostStep, 1);

  pManager->AddRestProcess(&thePionMinusAbsorption, ordDefault);

  // KaonPlus
  pManager = G4KaonPlus::KaonPlus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEKaonPlusModel = new G4LEKaonPlusInelastic();
  theHEKaonPlusModel = new G4HEKaonPlusInelastic();
  theKaonPlusInelastic.RegisterMe(theLEKaonPlusModel);
  theKaonPlusInelastic.RegisterMe(theHEKaonPlusModel);
  pManager->AddDiscreteProcess(&theKaonPlusInelastic);

  pManager->AddProcess(&theKaonPlusIonisation, ordInActive,2, 2);

  pManager->AddProcess(&theKaonPlusMult);
  pManager->SetProcessOrdering(&theKaonPlusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theKaonPlusMult, idxPostStep, 1);

  // KaonMinus
  pManager = G4KaonMinus::KaonMinus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEKaonMinusModel = new G4LEKaonMinusInelastic();
  theHEKaonMinusModel = new G4HEKaonMinusInelastic();
  theKaonMinusInelastic.RegisterMe(theLEKaonMinusModel);
  theKaonMinusInelastic.RegisterMe(theHEKaonMinusModel);
  pManager->AddDiscreteProcess(&theKaonMinusInelastic);

  pManager->AddProcess(&theKaonMinusIonisation, ordInActive,2, 2);

  pManager->AddProcess(&theKaonMinusMult);
  pManager->SetProcessOrdering(&theKaonMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theKaonMinusMult, idxPostStep, 1);

  pManager->AddRestProcess(&theKaonMinusAbsorption, ordDefault);

  // KaonZeroL
  pManager = G4KaonZeroLong::KaonZeroLong()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEKaonZeroLModel = new G4LEKaonZeroLInelastic();
  theHEKaonZeroLModel = new G4HEKaonZeroInelastic();
  theKaonZeroLInelastic.RegisterMe(theLEKaonZeroLModel);
  theKaonZeroLInelastic.RegisterMe(theHEKaonZeroLModel);
  pManager->AddDiscreteProcess(&theKaonZeroLInelastic);
 
  // KaonZeroS
  pManager = G4KaonZeroShort::KaonZeroShort()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEKaonZeroSModel = new G4LEKaonZeroSInelastic();
  theHEKaonZeroSModel = new G4HEKaonZeroInelastic();
  theKaonZeroSInelastic.RegisterMe(theLEKaonZeroSModel);
  theKaonZeroSInelastic.RegisterMe(theHEKaonZeroSModel);
  pManager->AddDiscreteProcess(&theKaonZeroSInelastic);

  // Proton
  pManager = G4Proton::Proton()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEProtonModel = new G4LEProtonInelastic();
  theHEProtonModel = new G4HEProtonInelastic();
  theProtonInelastic.RegisterMe(theLEProtonModel);
  theProtonInelastic.RegisterMe(theHEProtonModel);
  pManager->AddDiscreteProcess(&theProtonInelastic);

  pManager->AddProcess(&theProtonIonisation, ordInActive,2, 2);

  pManager->AddProcess(&theProtonMult);
  pManager->SetProcessOrdering(&theProtonMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theProtonMult, idxPostStep, 1);

  // anti-Proton
  pManager = G4AntiProton::AntiProton()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEAntiProtonModel = new G4LEAntiProtonInelastic();
  theHEAntiProtonModel = new G4HEAntiProtonInelastic();
  theAntiProtonInelastic.RegisterMe(theLEAntiProtonModel);
  theAntiProtonInelastic.RegisterMe(theHEAntiProtonModel);
  pManager->AddDiscreteProcess(&theAntiProtonInelastic);

  pManager->AddProcess(&theAntiProtonIonisation, ordInActive,2, 2);

  pManager->AddProcess(&theAntiProtonMult);
  pManager->SetProcessOrdering(&theAntiProtonMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theAntiProtonMult, idxPostStep, 1);

  pManager->AddRestProcess(&theAntiProtonAnnihilation);

  // Neutron
  pManager = G4Neutron::Neutron()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLENeutronModel = new G4LENeutronInelastic();
  theHENeutronModel = new G4HENeutronInelastic();
  theNeutronInelastic.RegisterMe(theLENeutronModel);
  theNeutronInelastic.RegisterMe(theHENeutronModel);
  pManager->AddDiscreteProcess(&theNeutronInelastic);
  
  //theNeutronFissionModel = new G4LFission();
  //theNeutronFission.RegisterMe(theNeutronFissionModel);
  //pManager->AddDiscreteProcess(&NeutronFission);

  //theNeutronCaptureModel = new G4LCapture();
  //theNeutronCapture.RegisterMe(theNeutronCaptureModel);
  //pManager->AddDiscreteProcess(&theNeutronCapture);

  // AntiNeutron
  pManager = G4AntiNeutron::AntiNeutron()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEAntiNeutronModel = new G4LEAntiNeutronInelastic();
  theHEAntiNeutronModel = new G4HEAntiNeutronInelastic();
  theAntiNeutronInelastic.RegisterMe(theLEAntiNeutronModel);
  theAntiNeutronInelastic.RegisterMe(theHEAntiNeutronModel);
  pManager->AddDiscreteProcess(&theAntiNeutronInelastic);
    
  pManager->AddRestProcess(&theAntiNeutronAnnihilation);

  // Lambda
  pManager = G4Lambda::Lambda()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLELambdaModel = new G4LELambdaInelastic();
  theHELambdaModel = new G4HELambdaInelastic();
  theLambdaInelastic.RegisterMe(theLELambdaModel);
  theLambdaInelastic.RegisterMe(theHELambdaModel);
  pManager->AddDiscreteProcess(&theLambdaInelastic);
  
  // AntiLambda
  pManager = G4AntiLambda::AntiLambda()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEAntiLambdaModel = new G4LEAntiLambdaInelastic();
  theHEAntiLambdaModel = new G4HEAntiLambdaInelastic();
  theAntiLambdaInelastic.RegisterMe(theLEAntiLambdaModel);
  theAntiLambdaInelastic.RegisterMe(theHEAntiLambdaModel);
  pManager->AddDiscreteProcess(&theAntiLambdaInelastic);
    
  // SigmaMinus
  pManager = G4SigmaMinus::SigmaMinus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLESigmaMinusModel = new G4LESigmaMinusInelastic();
  theHESigmaMinusModel = new G4HESigmaMinusInelastic();
  theSigmaMinusInelastic.RegisterMe(theLESigmaMinusModel);
  theSigmaMinusInelastic.RegisterMe(theHESigmaMinusModel);
  pManager->AddDiscreteProcess(&theSigmaMinusInelastic);

  pManager->AddProcess(&theSigmaMinusIonisation, ordInActive,2, 2);

  pManager->AddProcess(&theSigmaMinusMult);
  pManager->SetProcessOrdering(&theSigmaMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theSigmaMinusMult, idxPostStep, 1);

  // anti-SigmaMinus
  pManager = G4AntiSigmaMinus::AntiSigmaMinus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEAntiSigmaMinusModel = new G4LEAntiSigmaMinusInelastic();
  theHEAntiSigmaMinusModel = new G4HEAntiSigmaMinusInelastic();
  theAntiSigmaMinusInelastic.RegisterMe(theLEAntiSigmaMinusModel);
  theAntiSigmaMinusInelastic.RegisterMe(theHEAntiSigmaMinusModel);
  pManager->AddDiscreteProcess(&theAntiSigmaMinusInelastic);

  pManager->AddProcess(&theAntiSigmaMinusIonisation, ordInActive,2, 2);

  pManager->AddProcess(&theAntiSigmaMinusMult);
  pManager->SetProcessOrdering(&theAntiSigmaMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theAntiSigmaMinusMult, idxPostStep, 1);

  // SigmaPlus
  pManager = G4SigmaPlus::SigmaPlus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLESigmaPlusModel = new G4LESigmaPlusInelastic();
  theHESigmaPlusModel = new G4HESigmaPlusInelastic();
  theSigmaPlusInelastic.RegisterMe(theLESigmaPlusModel);
  theSigmaPlusInelastic.RegisterMe(theHESigmaPlusModel);
  pManager->AddDiscreteProcess(&theSigmaPlusInelastic);

  pManager->AddProcess(&theSigmaPlusIonisation, ordInActive,2, 2);

  pManager->AddProcess(&theSigmaPlusMult);
  pManager->SetProcessOrdering(&theSigmaPlusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theSigmaPlusMult, idxPostStep, 1);

  // anti-SigmaPlus
  pManager = G4AntiSigmaPlus::AntiSigmaPlus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEAntiSigmaPlusModel = new G4LEAntiSigmaPlusInelastic();
  theHEAntiSigmaPlusModel = new G4HEAntiSigmaPlusInelastic();
  theAntiSigmaPlusInelastic.RegisterMe(theLEAntiSigmaPlusModel);
  theAntiSigmaPlusInelastic.RegisterMe(theHEAntiSigmaPlusModel);
  pManager->AddDiscreteProcess(&theAntiSigmaPlusInelastic);

  pManager->AddProcess(&theAntiSigmaPlusIonisation, ordInActive,2, 2);

  pManager->AddProcess(&theAntiSigmaPlusMult);
  pManager->SetProcessOrdering(&theAntiSigmaPlusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theAntiSigmaPlusMult, idxPostStep, 1);

  // XiMinus
  pManager = G4XiMinus::XiMinus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEXiMinusModel = new G4LEXiMinusInelastic();
  theHEXiMinusModel = new G4HEXiMinusInelastic();
  theXiMinusInelastic.RegisterMe(theLEXiMinusModel);
  theXiMinusInelastic.RegisterMe(theHEXiMinusModel);
  pManager->AddDiscreteProcess(&theXiMinusInelastic);

  pManager->AddProcess(&theXiMinusIonisation, ordInActive,2, 2);

  pManager->AddProcess(&theXiMinusMult);
  pManager->SetProcessOrdering(&theXiMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theXiMinusMult, idxPostStep, 1);

  // anti-XiMinus
  pManager = G4AntiXiMinus::AntiXiMinus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEAntiXiMinusModel = new G4LEAntiXiMinusInelastic();
  theHEAntiXiMinusModel = new G4HEAntiXiMinusInelastic();
  theAntiXiMinusInelastic.RegisterMe(theLEAntiXiMinusModel);
  theAntiXiMinusInelastic.RegisterMe(theHEAntiXiMinusModel);
  pManager->AddDiscreteProcess(&theAntiXiMinusInelastic);

  pManager->AddProcess(&theAntiXiMinusIonisation, ordInActive,2, 2);

  pManager->AddProcess(&theAntiXiMinusMult);
  pManager->SetProcessOrdering(&theAntiXiMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theAntiXiMinusMult, idxPostStep, 1);

  // XiZero
  pManager = G4XiZero::XiZero()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEXiZeroModel = new G4LEXiZeroInelastic();
  theHEXiZeroModel = new G4HEXiZeroInelastic();
  theXiZeroInelastic.RegisterMe(theLEXiZeroModel);
  theXiZeroInelastic.RegisterMe(theHEXiZeroModel);
  pManager->AddDiscreteProcess(&theXiZeroInelastic);

  // anti-XiZero
  pManager = G4AntiXiZero::AntiXiZero()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEAntiXiZeroModel = new G4LEAntiXiZeroInelastic();
  theHEAntiXiZeroModel = new G4HEAntiXiZeroInelastic();
  theAntiXiZeroInelastic.RegisterMe(theLEAntiXiZeroModel);
  theAntiXiZeroInelastic.RegisterMe(theHEAntiXiZeroModel);
  pManager->AddDiscreteProcess(&theAntiXiZeroInelastic);

  // OmegaMinus
  pManager = G4OmegaMinus::OmegaMinus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEOmegaMinusModel = new G4LEOmegaMinusInelastic();
  theHEOmegaMinusModel = new G4HEOmegaMinusInelastic();
  theOmegaMinusInelastic.RegisterMe(theLEOmegaMinusModel);
  theOmegaMinusInelastic.RegisterMe(theHEOmegaMinusModel);
  pManager->AddDiscreteProcess(&theOmegaMinusInelastic);

  pManager->AddProcess(&theOmegaMinusIonisation, ordInActive,2, 2);

  pManager->AddProcess(&theOmegaMinusMult);
  pManager->SetProcessOrdering(&theOmegaMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theOmegaMinusMult, idxPostStep, 1);

  // anti-OmegaMinus
  pManager = G4AntiOmegaMinus::AntiOmegaMinus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEAntiOmegaMinusModel = new G4LEAntiOmegaMinusInelastic();
  theHEAntiOmegaMinusModel = new G4HEAntiOmegaMinusInelastic();
  theAntiOmegaMinusInelastic.RegisterMe(theLEAntiOmegaMinusModel);
  theAntiOmegaMinusInelastic.RegisterMe(theHEAntiOmegaMinusModel);
  pManager->AddDiscreteProcess(&theAntiOmegaMinusInelastic);

  pManager->AddProcess(&theAntiOmegaMinusIonisation, ordInActive,2, 2);

  pManager->AddProcess(&theAntiOmegaMinusMult);
  pManager->SetProcessOrdering(&theAntiOmegaMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theAntiOmegaMinusMult, idxPostStep, 1);

}





