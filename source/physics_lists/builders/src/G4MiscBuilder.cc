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
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include "G4MiscBuilder.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

#include "G4ChipsHyperonInelasticXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4ComponentAntiNuclNuclearXS.hh"


G4MiscBuilder::G4MiscBuilder(): 
  theAntiProtonInelastic(0), theLEAntiProtonModel(0), 
  theHEAntiProtonModel(0),
  theAntiNeutronInelastic(0), theLEAntiNeutronModel(0),
  theHEAntiNeutronModel(0),
  theLambdaInelastic(0), theLELambdaModel(0), theHELambdaModel(0),
  theAntiLambdaInelastic(0), theLEAntiLambdaModel(0), theHEAntiLambdaModel(0),
  theSigmaMinusInelastic(0), theLESigmaMinusModel(0), theHESigmaMinusModel(0),
  theAntiSigmaMinusInelastic(0), theLEAntiSigmaMinusModel(0), theHEAntiSigmaMinusModel(0),
  theSigmaPlusInelastic(0), theLESigmaPlusModel(0), theHESigmaPlusModel(0),
  theAntiSigmaPlusInelastic(0), theLEAntiSigmaPlusModel(0), theHEAntiSigmaPlusModel(0),
  theXiZeroInelastic(0), theLEXiZeroModel(0), theHEXiZeroModel(0),
  theAntiXiZeroInelastic(0), theLEAntiXiZeroModel(0), theHEAntiXiZeroModel(0),
  theXiMinusInelastic(0), theLEXiMinusModel(0), theHEXiMinusModel(0),
  theAntiXiMinusInelastic(0), theLEAntiXiMinusModel(0), theHEAntiXiMinusModel(0),
  theOmegaMinusInelastic(0), theLEOmegaMinusModel(0), theHEOmegaMinusModel(0),
  theAntiOmegaMinusInelastic(0), theLEAntiOmegaMinusModel(0), theHEAntiOmegaMinusModel(0),
  wasActivated(false)

{}


G4MiscBuilder::~G4MiscBuilder()
{
}

void G4MiscBuilder::Build()
{
  G4ProcessManager * aProcMan = 0;
  wasActivated = true;
  
  theAntiNucleonData = new G4CrossSectionInelastic(new G4ComponentAntiNuclNuclearXS());
  theChipsInelastic = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsHyperonInelasticXS::Default_Name());
    
  // anti-Proton
  theAntiProtonInelastic = new G4AntiProtonInelasticProcess();
  aProcMan = G4AntiProton::AntiProton()->GetProcessManager();
  theLEAntiProtonModel = new G4LEAntiProtonInelastic();
  theHEAntiProtonModel = new G4HEAntiProtonInelastic();
  theHEAntiProtonModel->SetMaxEnergy(100*TeV);
  theAntiProtonInelastic->RegisterMe(theLEAntiProtonModel);
  theAntiProtonInelastic->RegisterMe(theHEAntiProtonModel);
  aProcMan->AddDiscreteProcess(theAntiProtonInelastic);
  theAntiProtonInelastic->AddDataSet(theAntiNucleonData);

  // AntiNeutron
  theAntiNeutronInelastic = new G4AntiNeutronInelasticProcess();
  aProcMan = G4AntiNeutron::AntiNeutron()->GetProcessManager();
  theLEAntiNeutronModel = new G4LEAntiNeutronInelastic();
  theHEAntiNeutronModel = new G4HEAntiNeutronInelastic();
  theHEAntiNeutronModel->SetMaxEnergy(100*TeV);
  theAntiNeutronInelastic->RegisterMe(theLEAntiNeutronModel);
  theAntiNeutronInelastic->RegisterMe(theHEAntiNeutronModel);
  aProcMan->AddDiscreteProcess(theAntiNeutronInelastic);
  theAntiNeutronInelastic->AddDataSet(theAntiNucleonData);
    
  // Lambda
  theLambdaInelastic = new G4LambdaInelasticProcess();
  aProcMan = G4Lambda::Lambda()->GetProcessManager();
  theLELambdaModel = new G4LELambdaInelastic();
  theHELambdaModel = new G4HELambdaInelastic();
  theHELambdaModel->SetMaxEnergy(100*TeV);
  theLambdaInelastic->RegisterMe(theLELambdaModel);
  theLambdaInelastic->RegisterMe(theHELambdaModel);
  aProcMan->AddDiscreteProcess(theLambdaInelastic);
  theLambdaInelastic->AddDataSet(theChipsInelastic);
    
  // AntiLambda
  theAntiLambdaInelastic = new G4AntiLambdaInelasticProcess();
  aProcMan = G4AntiLambda::AntiLambda()->GetProcessManager();
  theLEAntiLambdaModel = new G4LEAntiLambdaInelastic();
  theHEAntiLambdaModel = new G4HEAntiLambdaInelastic();
  theHEAntiLambdaModel->SetMaxEnergy(100*TeV);
  theAntiLambdaInelastic->RegisterMe(theLEAntiLambdaModel);
  theAntiLambdaInelastic->RegisterMe(theHEAntiLambdaModel);
  aProcMan->AddDiscreteProcess(theAntiLambdaInelastic);
  theAntiLambdaInelastic->AddDataSet(theChipsInelastic);
    
  // SigmaMinus
  theSigmaMinusInelastic = new G4SigmaMinusInelasticProcess();
  aProcMan = G4SigmaMinus::SigmaMinus()->GetProcessManager();
  theLESigmaMinusModel = new G4LESigmaMinusInelastic();
  theHESigmaMinusModel = new G4HESigmaMinusInelastic();
  theHESigmaMinusModel->SetMaxEnergy(100*TeV);
  theSigmaMinusInelastic->RegisterMe(theLESigmaMinusModel);
  theSigmaMinusInelastic->RegisterMe(theHESigmaMinusModel);
  aProcMan->AddDiscreteProcess(theSigmaMinusInelastic);
  theSigmaMinusInelastic->AddDataSet(theChipsInelastic);
    
  // anti-SigmaMinus
  theAntiSigmaMinusInelastic = new G4AntiSigmaMinusInelasticProcess();
  aProcMan = G4AntiSigmaMinus::AntiSigmaMinus()->GetProcessManager();
  theLEAntiSigmaMinusModel = new G4LEAntiSigmaMinusInelastic();
  theHEAntiSigmaMinusModel = new G4HEAntiSigmaMinusInelastic();
  theHEAntiSigmaMinusModel->SetMaxEnergy(100*TeV);
  theAntiSigmaMinusInelastic->RegisterMe(theLEAntiSigmaMinusModel);
  theAntiSigmaMinusInelastic->RegisterMe(theHEAntiSigmaMinusModel);
  aProcMan->AddDiscreteProcess(theAntiSigmaMinusInelastic);
  theAntiSigmaMinusInelastic->AddDataSet(theChipsInelastic);

  // SigmaPlus
  theSigmaPlusInelastic = new G4SigmaPlusInelasticProcess();
  aProcMan = G4SigmaPlus::SigmaPlus()->GetProcessManager();
  theLESigmaPlusModel = new G4LESigmaPlusInelastic();
  theHESigmaPlusModel = new G4HESigmaPlusInelastic();
  theHESigmaPlusModel->SetMaxEnergy(100*TeV);
  theSigmaPlusInelastic->RegisterMe(theLESigmaPlusModel);
  theSigmaPlusInelastic->RegisterMe(theHESigmaPlusModel);
  aProcMan->AddDiscreteProcess(theSigmaPlusInelastic);
  theSigmaPlusInelastic->AddDataSet(theChipsInelastic);
    
  // anti-SigmaPlus
  theAntiSigmaPlusInelastic = new G4AntiSigmaPlusInelasticProcess();
  aProcMan = G4AntiSigmaPlus::AntiSigmaPlus()->GetProcessManager();
  theLEAntiSigmaPlusModel = new G4LEAntiSigmaPlusInelastic();
  theHEAntiSigmaPlusModel = new G4HEAntiSigmaPlusInelastic();
  theHEAntiSigmaPlusModel->SetMaxEnergy(100*TeV);
  theAntiSigmaPlusInelastic->RegisterMe(theLEAntiSigmaPlusModel);
  theAntiSigmaPlusInelastic->RegisterMe(theHEAntiSigmaPlusModel);
  aProcMan->AddDiscreteProcess(theAntiSigmaPlusInelastic);
  theAntiSigmaPlusInelastic->AddDataSet(theChipsInelastic);

  // XiMinus
  theXiMinusInelastic = new G4XiMinusInelasticProcess();
  aProcMan = G4XiMinus::XiMinus()->GetProcessManager();
  theLEXiMinusModel = new G4LEXiMinusInelastic();
  theHEXiMinusModel = new G4HEXiMinusInelastic();
  theHEXiMinusModel->SetMaxEnergy(100*TeV);
  theXiMinusInelastic->RegisterMe(theLEXiMinusModel);
  theXiMinusInelastic->RegisterMe(theHEXiMinusModel);
  aProcMan->AddDiscreteProcess(theXiMinusInelastic);
  theXiMinusInelastic->AddDataSet(theChipsInelastic);

  // anti-XiMinus
  theAntiXiMinusInelastic = new G4AntiXiMinusInelasticProcess();
  aProcMan = G4AntiXiMinus::AntiXiMinus()->GetProcessManager();
  theLEAntiXiMinusModel = new G4LEAntiXiMinusInelastic();
  theHEAntiXiMinusModel = new G4HEAntiXiMinusInelastic();
  theHEAntiXiMinusModel->SetMaxEnergy(100*TeV);
  theAntiXiMinusInelastic->RegisterMe(theLEAntiXiMinusModel);
  theAntiXiMinusInelastic->RegisterMe(theHEAntiXiMinusModel);
  aProcMan->AddDiscreteProcess(theAntiXiMinusInelastic);
  theAntiXiMinusInelastic->AddDataSet(theChipsInelastic);

  // XiZero
  theXiZeroInelastic = new G4XiZeroInelasticProcess();
  aProcMan = G4XiZero::XiZero()->GetProcessManager();
  theLEXiZeroModel = new G4LEXiZeroInelastic();
  theHEXiZeroModel = new G4HEXiZeroInelastic();
  theHEXiZeroModel->SetMaxEnergy(100*TeV);
  theXiZeroInelastic->RegisterMe(theLEXiZeroModel);
  theXiZeroInelastic->RegisterMe(theHEXiZeroModel);
  aProcMan->AddDiscreteProcess(theXiZeroInelastic);
  theXiZeroInelastic->AddDataSet(theChipsInelastic);
    
  // anti-XiZero
  theAntiXiZeroInelastic = new G4AntiXiZeroInelasticProcess();
  aProcMan = G4AntiXiZero::AntiXiZero()->GetProcessManager();
  theLEAntiXiZeroModel = new G4LEAntiXiZeroInelastic();
  theHEAntiXiZeroModel = new G4HEAntiXiZeroInelastic();
  theHEAntiXiZeroModel->SetMaxEnergy(100*TeV);
  theAntiXiZeroInelastic->RegisterMe(theLEAntiXiZeroModel);
  theAntiXiZeroInelastic->RegisterMe(theHEAntiXiZeroModel);
  aProcMan->AddDiscreteProcess(theAntiXiZeroInelastic);
  theAntiXiZeroInelastic->AddDataSet(theChipsInelastic);

  // OmegaMinus
  theOmegaMinusInelastic = new G4OmegaMinusInelasticProcess();
  aProcMan = G4OmegaMinus::OmegaMinus()->GetProcessManager();
  theLEOmegaMinusModel = new G4LEOmegaMinusInelastic();
  theHEOmegaMinusModel = new G4HEOmegaMinusInelastic();
  theHEOmegaMinusModel->SetMaxEnergy(100*TeV);
  theOmegaMinusInelastic->RegisterMe(theLEOmegaMinusModel);
  theOmegaMinusInelastic->RegisterMe(theHEOmegaMinusModel);
  aProcMan->AddDiscreteProcess(theOmegaMinusInelastic);
  theOmegaMinusInelastic->AddDataSet(theChipsInelastic);

  // anti-OmegaMinus
  theAntiOmegaMinusInelastic = new G4AntiOmegaMinusInelasticProcess();
  aProcMan = G4AntiOmegaMinus::AntiOmegaMinus()->GetProcessManager();
  theLEAntiOmegaMinusModel = new G4LEAntiOmegaMinusInelastic();
  theHEAntiOmegaMinusModel = new G4HEAntiOmegaMinusInelastic();
  theHEAntiOmegaMinusModel->SetMaxEnergy(100*TeV);
  theAntiOmegaMinusInelastic->RegisterMe(theLEAntiOmegaMinusModel);
  theAntiOmegaMinusInelastic->RegisterMe(theHEAntiOmegaMinusModel);
  aProcMan->AddDiscreteProcess(theAntiOmegaMinusInelastic);
  theAntiOmegaMinusInelastic->AddDataSet(theChipsInelastic);
}

