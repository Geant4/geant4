#include "G4MiscLHEPBuilder.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4MiscLHEPBuilder::G4MiscLHEPBuilder() {}
G4MiscLHEPBuilder::~G4MiscLHEPBuilder() {}

void G4MiscLHEPBuilder::Build()
{
  G4ProcessManager * aProcMan = 0;
  theElasticModel = new G4LElastic;
  
  // anti-Proton
  aProcMan = G4AntiProton::AntiProton()->GetProcessManager();
  theLEAntiProtonModel = new G4LEAntiProtonInelastic();
  theHEAntiProtonModel = new G4HEAntiProtonInelastic();
  theAntiProtonInelastic.RegisterMe(theLEAntiProtonModel);
  theAntiProtonInelastic.RegisterMe(theHEAntiProtonModel);
  aProcMan->AddDiscreteProcess(&theAntiProtonInelastic);
  theAntiProtonElasticProcess.RegisterMe(theElasticModel);
  aProcMan->AddDiscreteProcess(&theAntiProtonElasticProcess);

  // AntiNeutron
  aProcMan = G4AntiNeutron::AntiNeutron()->GetProcessManager();
  theLEAntiNeutronModel = new G4LEAntiNeutronInelastic();
  theHEAntiNeutronModel = new G4HEAntiNeutronInelastic();
  theAntiNeutronInelastic.RegisterMe(theLEAntiNeutronModel);
  theAntiNeutronInelastic.RegisterMe(theHEAntiNeutronModel);
  aProcMan->AddDiscreteProcess(&theAntiNeutronInelastic);
  theAntiNeutronElasticProcess.RegisterMe(theElasticModel);
  aProcMan->AddDiscreteProcess(&theAntiNeutronElasticProcess);

  // Lambda
  aProcMan = G4Lambda::Lambda()->GetProcessManager();
  theLELambdaModel = new G4LELambdaInelastic();
  theHELambdaModel = new G4HELambdaInelastic();
  theLambdaInelastic.RegisterMe(theLELambdaModel);
  theLambdaInelastic.RegisterMe(theHELambdaModel);
  aProcMan->AddDiscreteProcess(&theLambdaInelastic);
  theLambdaElasticProcess.RegisterMe(theElasticModel);
  aProcMan->AddDiscreteProcess(&theLambdaElasticProcess);
  
  // AntiLambda
  aProcMan = G4AntiLambda::AntiLambda()->GetProcessManager();
  theLEAntiLambdaModel = new G4LEAntiLambdaInelastic();
  theHEAntiLambdaModel = new G4HEAntiLambdaInelastic();
  theAntiLambdaInelastic.RegisterMe(theLEAntiLambdaModel);
  theAntiLambdaInelastic.RegisterMe(theHEAntiLambdaModel);
  aProcMan->AddDiscreteProcess(&theAntiLambdaInelastic);
  theAntiLambdaElasticProcess.RegisterMe(theElasticModel);
  aProcMan->AddDiscreteProcess(&theAntiLambdaElasticProcess);
    
  // SigmaMinus
  aProcMan = G4SigmaMinus::SigmaMinus()->GetProcessManager();
  theLESigmaMinusModel = new G4LESigmaMinusInelastic();
  theHESigmaMinusModel = new G4HESigmaMinusInelastic();
  theSigmaMinusInelastic.RegisterMe(theLESigmaMinusModel);
  theSigmaMinusInelastic.RegisterMe(theHESigmaMinusModel);
  aProcMan->AddDiscreteProcess(&theSigmaMinusInelastic);
  theSigmaMinusElasticProcess.RegisterMe(theElasticModel);
  aProcMan->AddDiscreteProcess(&theSigmaMinusElasticProcess);

  // anti-SigmaMinus
  aProcMan = G4AntiSigmaMinus::AntiSigmaMinus()->GetProcessManager();
  theLEAntiSigmaMinusModel = new G4LEAntiSigmaMinusInelastic();
  theHEAntiSigmaMinusModel = new G4HEAntiSigmaMinusInelastic();
  theAntiSigmaMinusInelastic.RegisterMe(theLEAntiSigmaMinusModel);
  theAntiSigmaMinusInelastic.RegisterMe(theHEAntiSigmaMinusModel);
  aProcMan->AddDiscreteProcess(&theAntiSigmaMinusInelastic);
  theAntiSigmaMinusElasticProcess.RegisterMe(theElasticModel);
  aProcMan->AddDiscreteProcess(&theAntiSigmaMinusElasticProcess);

  // SigmaPlus
  aProcMan = G4SigmaPlus::SigmaPlus()->GetProcessManager();
  theLESigmaPlusModel = new G4LESigmaPlusInelastic();
  theHESigmaPlusModel = new G4HESigmaPlusInelastic();
  theSigmaPlusInelastic.RegisterMe(theLESigmaPlusModel);
  theSigmaPlusInelastic.RegisterMe(theHESigmaPlusModel);
  aProcMan->AddDiscreteProcess(&theSigmaPlusInelastic);
  theSigmaPlusElasticProcess.RegisterMe(theElasticModel);
  aProcMan->AddDiscreteProcess(&theSigmaPlusElasticProcess);

  // anti-SigmaPlus
  aProcMan = G4AntiSigmaPlus::AntiSigmaPlus()->GetProcessManager();
  theLEAntiSigmaPlusModel = new G4LEAntiSigmaPlusInelastic();
  theHEAntiSigmaPlusModel = new G4HEAntiSigmaPlusInelastic();
  theAntiSigmaPlusInelastic.RegisterMe(theLEAntiSigmaPlusModel);
  theAntiSigmaPlusInelastic.RegisterMe(theHEAntiSigmaPlusModel);
  aProcMan->AddDiscreteProcess(&theAntiSigmaPlusInelastic);
  theAntiSigmaPlusElasticProcess.RegisterMe(theElasticModel);
  aProcMan->AddDiscreteProcess(&theAntiSigmaPlusElasticProcess);

  // XiMinus
  aProcMan = G4XiMinus::XiMinus()->GetProcessManager();
  theLEXiMinusModel = new G4LEXiMinusInelastic();
  theHEXiMinusModel = new G4HEXiMinusInelastic();
  theXiMinusInelastic.RegisterMe(theLEXiMinusModel);
  theXiMinusInelastic.RegisterMe(theHEXiMinusModel);
  aProcMan->AddDiscreteProcess(&theXiMinusInelastic);
  theXiMinusElasticProcess.RegisterMe(theElasticModel);
  aProcMan->AddDiscreteProcess(&theXiMinusElasticProcess);

  // anti-XiMinus
  aProcMan = G4AntiXiMinus::AntiXiMinus()->GetProcessManager();
  theLEAntiXiMinusModel = new G4LEAntiXiMinusInelastic();
  theHEAntiXiMinusModel = new G4HEAntiXiMinusInelastic();
  theAntiXiMinusInelastic.RegisterMe(theLEAntiXiMinusModel);
  theAntiXiMinusInelastic.RegisterMe(theHEAntiXiMinusModel);
  aProcMan->AddDiscreteProcess(&theAntiXiMinusInelastic);
  theAntiXiMinusElasticProcess.RegisterMe(theElasticModel);
  aProcMan->AddDiscreteProcess(&theAntiXiMinusElasticProcess);

  // XiZero
  aProcMan = G4XiZero::XiZero()->GetProcessManager();
  theLEXiZeroModel = new G4LEXiZeroInelastic();
  theHEXiZeroModel = new G4HEXiZeroInelastic();
  theXiZeroInelastic.RegisterMe(theLEXiZeroModel);
  theXiZeroInelastic.RegisterMe(theHEXiZeroModel);
  aProcMan->AddDiscreteProcess(&theXiZeroInelastic);
  theXiZeroElasticProcess.RegisterMe(theElasticModel);
  aProcMan->AddDiscreteProcess(&theXiZeroElasticProcess);

  // anti-XiZero
  aProcMan = G4AntiXiZero::AntiXiZero()->GetProcessManager();
  theLEAntiXiZeroModel = new G4LEAntiXiZeroInelastic();
  theHEAntiXiZeroModel = new G4HEAntiXiZeroInelastic();
  theAntiXiZeroInelastic.RegisterMe(theLEAntiXiZeroModel);
  theAntiXiZeroInelastic.RegisterMe(theHEAntiXiZeroModel);
  aProcMan->AddDiscreteProcess(&theAntiXiZeroInelastic);
  theAntiXiZeroElasticProcess.RegisterMe(theElasticModel);
  aProcMan->AddDiscreteProcess(&theAntiXiZeroElasticProcess);

  // OmegaMinus
  aProcMan = G4OmegaMinus::OmegaMinus()->GetProcessManager();
  theLEOmegaMinusModel = new G4LEOmegaMinusInelastic();
  theHEOmegaMinusModel = new G4HEOmegaMinusInelastic();
  theOmegaMinusInelastic.RegisterMe(theLEOmegaMinusModel);
  theOmegaMinusInelastic.RegisterMe(theHEOmegaMinusModel);
  aProcMan->AddDiscreteProcess(&theOmegaMinusInelastic);
  theOmegaMinusElasticProcess.RegisterMe(theElasticModel);
  aProcMan->AddDiscreteProcess(&theOmegaMinusElasticProcess);

  // anti-OmegaMinus
  aProcMan = G4AntiOmegaMinus::AntiOmegaMinus()->GetProcessManager();
  theLEAntiOmegaMinusModel = new G4LEAntiOmegaMinusInelastic();
  theHEAntiOmegaMinusModel = new G4HEAntiOmegaMinusInelastic();
  theAntiOmegaMinusInelastic.RegisterMe(theLEAntiOmegaMinusModel);
  theAntiOmegaMinusInelastic.RegisterMe(theHEAntiOmegaMinusModel);
  aProcMan->AddDiscreteProcess(&theAntiOmegaMinusInelastic);
  theAntiOmegaMinusElasticProcess.RegisterMe(theElasticModel);
  aProcMan->AddDiscreteProcess(&theAntiOmegaMinusElasticProcess);
}

// 2002 by J.P. Wellisch
