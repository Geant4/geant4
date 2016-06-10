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

#include "ExExChHyperonFTFPBuilder.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4CrossSectionDataSetRegistry.hh"

// Wrapper
#include "XWrapperDiscreteProcess.hh"

ExExChHyperonFTFPBuilder::ExExChHyperonFTFPBuilder():
theLambdaInelastic(0),
theAntiLambdaInelastic(0),
theSigmaMinusInelastic(0),
theAntiSigmaMinusInelastic(0),
theSigmaPlusInelastic(0),
theAntiSigmaPlusInelastic(0),
theXiZeroInelastic(0),
theAntiXiZeroInelastic(0),
theXiMinusInelastic(0),
theAntiXiMinusInelastic(0),
theOmegaMinusInelastic(0),
theAntiOmegaMinusInelastic(0),
wasActivated(false)
{
    
    // Hyperon : Bertini at low energies, then FTFP
    
    HyperonFTFP = new G4TheoFSGenerator("FTFP");
    
    HyperonFTFP->SetMinEnergy( 2.*GeV );
    HyperonFTFP->SetMaxEnergy( 100.*TeV );
    
    theStringModel = new G4FTFModel;
    theLund = new G4LundStringFragmentation;
    theStringDecay = new G4ExcitedStringDecay(theLund);
    theStringModel->SetFragmentationModel(theStringDecay);
    
    theCascade = new G4GeneratorPrecompoundInterface;
    
    theHandler =new G4ExcitationHandler;
    thePreEquilib = new G4PreCompoundModel(theHandler);
    theCascade->SetDeExcitation(thePreEquilib);
    
    HyperonFTFP->SetTransport(theCascade);
    HyperonFTFP->SetHighEnergyGenerator(theStringModel);
    
    theBertini = new G4CascadeInterface;
    theBertini->SetMinEnergy( 0.*GeV );
    theBertini->SetMaxEnergy( 6.*GeV );
    
    // AntiHyperons: Use FTFP for full energy range, starting at 0.
    
    AntiHyperonFTFP = new G4TheoFSGenerator("FTFP");
    AntiHyperonFTFP->SetMinEnergy( 0.*GeV );
    AntiHyperonFTFP->SetMaxEnergy( 100.*TeV );
    AntiHyperonFTFP->SetTransport(theCascade);
    AntiHyperonFTFP->SetHighEnergyGenerator(theStringModel);
    
    // use CHIPS cross sections
    theCHIPSInelastic = G4CrossSectionDataSetRegistry::
        Instance()->GetCrossSectionDataSet(
        G4ChipsHyperonInelasticXS::Default_Name());
}


ExExChHyperonFTFPBuilder::~ExExChHyperonFTFPBuilder()
{
    delete HyperonFTFP;
    delete theStringModel;
    delete theStringDecay;
    delete theCascade;
    delete thePreEquilib;
    //  delete theHandler;
    delete theBertini;
    delete AntiHyperonFTFP;
    
    if (wasActivated) {
        delete theLambdaInelastic;
        delete theAntiLambdaInelastic;
        delete theSigmaMinusInelastic;
        delete theAntiSigmaMinusInelastic;
        delete theSigmaPlusInelastic;
        delete theAntiSigmaPlusInelastic;
        delete theXiMinusInelastic;
        delete theAntiXiMinusInelastic;
        delete theXiZeroInelastic;
        delete theAntiXiZeroInelastic;
        delete theOmegaMinusInelastic;
        delete theAntiOmegaMinusInelastic;
    }
}

void ExExChHyperonFTFPBuilder::Build()
{
    G4ProcessManager * aProcMan = 0;
    wasActivated = true;
    
    // Lambda
    theLambdaInelastic = new G4LambdaInelasticProcess();
    theLambdaInelastic->RegisterMe(theBertini);
    theLambdaInelastic->RegisterMe(HyperonFTFP);
    theLambdaInelastic->AddDataSet(theCHIPSInelastic);
    aProcMan = G4Lambda::Lambda()->GetProcessManager();
    aProcMan->AddDiscreteProcess(theLambdaInelastic);
    
    // AntiLambda
    theAntiLambdaInelastic = new G4AntiLambdaInelasticProcess();
    theAntiLambdaInelastic->RegisterMe(AntiHyperonFTFP);
    theAntiLambdaInelastic->AddDataSet(theCHIPSInelastic);
    
    aProcMan = G4AntiLambda::AntiLambda()->GetProcessManager();
    aProcMan->AddDiscreteProcess(theAntiLambdaInelastic);
    
    // SigmaMinus
    theSigmaMinusInelastic = new G4SigmaMinusInelasticProcess();
    theSigmaMinusInelastic->RegisterMe(theBertini);
    theSigmaMinusInelastic->RegisterMe(HyperonFTFP);
    theSigmaMinusInelastic->AddDataSet(theCHIPSInelastic);
    
    aProcMan = G4SigmaMinus::SigmaMinus()->GetProcessManager();
    XWrapperDiscreteProcess* theSigmaMinusInelastic_w =
        new XWrapperDiscreteProcess();
    theSigmaMinusInelastic_w->RegisterProcess(theSigmaMinusInelastic,1);
    aProcMan->AddDiscreteProcess(theSigmaMinusInelastic_w);
    
    // anti-SigmaMinus
    theAntiSigmaMinusInelastic = new G4AntiSigmaMinusInelasticProcess();
    theAntiSigmaMinusInelastic->RegisterMe(AntiHyperonFTFP);
    theAntiSigmaMinusInelastic->AddDataSet(theCHIPSInelastic);
    
    aProcMan = G4AntiSigmaMinus::AntiSigmaMinus()->GetProcessManager();
    XWrapperDiscreteProcess* theAntiSigmaMinusInelastic_w =
        new XWrapperDiscreteProcess();
    theAntiSigmaMinusInelastic_w->RegisterProcess(theAntiSigmaMinusInelastic,1);
    aProcMan->AddDiscreteProcess(theAntiSigmaMinusInelastic_w);
    
    // SigmaPlus
    theSigmaPlusInelastic = new G4SigmaPlusInelasticProcess();
    theSigmaPlusInelastic->RegisterMe(theBertini);
    theSigmaPlusInelastic->RegisterMe(HyperonFTFP);
    theSigmaPlusInelastic->AddDataSet(theCHIPSInelastic);
    
    aProcMan = G4SigmaPlus::SigmaPlus()->GetProcessManager();
    XWrapperDiscreteProcess* theSigmaPlusInelastic_w =
        new XWrapperDiscreteProcess();
    theSigmaPlusInelastic_w->RegisterProcess(theSigmaPlusInelastic,1);
    aProcMan->AddDiscreteProcess(theSigmaPlusInelastic_w);
    
    // anti-SigmaPlus
    theAntiSigmaPlusInelastic = new G4AntiSigmaPlusInelasticProcess();
    theAntiSigmaPlusInelastic->RegisterMe(AntiHyperonFTFP);
    theAntiSigmaPlusInelastic->AddDataSet(theCHIPSInelastic);
    
    aProcMan = G4AntiSigmaPlus::AntiSigmaPlus()->GetProcessManager();
    XWrapperDiscreteProcess* theAntiSigmaPlusInelastic_w =
        new XWrapperDiscreteProcess();
    theAntiSigmaPlusInelastic_w->RegisterProcess(theAntiSigmaPlusInelastic,1);
    aProcMan->AddDiscreteProcess(theAntiSigmaPlusInelastic_w);
    
    // XiMinus
    theXiMinusInelastic = new G4XiMinusInelasticProcess();
    theXiMinusInelastic->RegisterMe(theBertini);
    theXiMinusInelastic->RegisterMe(HyperonFTFP);
    theXiMinusInelastic->AddDataSet(theCHIPSInelastic);
    
    aProcMan = G4XiMinus::XiMinus()->GetProcessManager();
    XWrapperDiscreteProcess* theXiMinusInelastic_w =
        new XWrapperDiscreteProcess();
    theXiMinusInelastic_w->RegisterProcess(theXiMinusInelastic,1);
    aProcMan->AddDiscreteProcess(theXiMinusInelastic_w);
    
    // anti-XiMinus
    theAntiXiMinusInelastic = new G4AntiXiMinusInelasticProcess();
    theAntiXiMinusInelastic->RegisterMe(AntiHyperonFTFP);
    theAntiXiMinusInelastic->AddDataSet(theCHIPSInelastic);
    
    aProcMan = G4AntiXiMinus::AntiXiMinus()->GetProcessManager();
    XWrapperDiscreteProcess* theAntiXiMinusInelastic_w =
        new XWrapperDiscreteProcess();
    theAntiXiMinusInelastic_w->RegisterProcess(theAntiXiMinusInelastic,1);
    aProcMan->AddDiscreteProcess(theAntiXiMinusInelastic_w);
    
    // XiZero
    theXiZeroInelastic = new G4XiZeroInelasticProcess();
    theXiZeroInelastic->RegisterMe(theBertini);
    theXiZeroInelastic->RegisterMe(HyperonFTFP);
    theXiZeroInelastic->AddDataSet(theCHIPSInelastic);
    
    aProcMan = G4XiZero::XiZero()->GetProcessManager();
    aProcMan->AddDiscreteProcess(theXiZeroInelastic);
    
    // anti-XiZero
    theAntiXiZeroInelastic = new G4AntiXiZeroInelasticProcess();
    theAntiXiZeroInelastic->RegisterMe(AntiHyperonFTFP);
    theAntiXiZeroInelastic->AddDataSet(theCHIPSInelastic);
    
    aProcMan = G4AntiXiZero::AntiXiZero()->GetProcessManager();
    aProcMan->AddDiscreteProcess(theAntiXiZeroInelastic);
    
    // OmegaMinus
    theOmegaMinusInelastic = new G4OmegaMinusInelasticProcess();
    theOmegaMinusInelastic->RegisterMe(theBertini);
    theOmegaMinusInelastic->RegisterMe(HyperonFTFP);
    theOmegaMinusInelastic->AddDataSet(theCHIPSInelastic);
    
    aProcMan = G4OmegaMinus::OmegaMinus()->GetProcessManager();
    XWrapperDiscreteProcess* theOmegaMinusInelastic_w =
        new XWrapperDiscreteProcess();
    theOmegaMinusInelastic_w->RegisterProcess(theOmegaMinusInelastic,1);
    aProcMan->AddDiscreteProcess(theOmegaMinusInelastic_w);
    
    // anti-OmegaMinus
    theAntiOmegaMinusInelastic = new G4AntiOmegaMinusInelasticProcess();
    theAntiOmegaMinusInelastic->RegisterMe(AntiHyperonFTFP);
    theAntiOmegaMinusInelastic->AddDataSet(theCHIPSInelastic);
    
    aProcMan = G4AntiOmegaMinus::AntiOmegaMinus()->GetProcessManager();
    XWrapperDiscreteProcess* theAntiOmegaMinusInelastic_w =
        new XWrapperDiscreteProcess();
    theAntiOmegaMinusInelastic_w->RegisterProcess(theAntiOmegaMinusInelastic,1);
    aProcMan->AddDiscreteProcess(theAntiOmegaMinusInelastic_w);
}


