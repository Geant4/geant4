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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: LowEHadronNuclear.cc,v 1.1 2001-09-21 17:38:47 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Johannes Peter Wellisch, 22.Apr 1997: full test-suite coded.    
    
#include "HadronNuclear/batch_monitor/ParticleInfo.h"

#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4Timer.hh"
 
#include "G4Material.hh"

#include "G4BosonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
 
#include "G4ProcessManager.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4AntiLambdaInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4AntiOmegaMinusInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4AntiSigmaMinusInelasticProcess.hh"
#include "G4AntiSigmaPlusInelasticProcess.hh"
#include "G4AntiXiMinusInelasticProcess.hh"
#include "G4AntiXiZeroInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4LambdaInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4OmegaMinusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4SigmaMinusInelasticProcess.hh"
#include "G4SigmaPlusInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4XiMinusInelasticProcess.hh"
#include "G4XiZeroInelasticProcess.hh"
 

#include "G4DynamicParticle.hh"
#include "G4Proton.hh"
#include "G4HardronNuclearChips.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4PreCompoundModel.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4GRSVolume.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"


 int main()
  {
    
    RanecuEngine theEngine;
    HepRandom::setTheEngine( &theEngine );
    theEngine.showStatus();
    G4cout << G4endl;
    
    G4String name, symbol;
    G4double a, iz, z, density;
    G4int nEl;
    
    G4Material *thePb = new G4Material(name="Lead", density=11.35*g/cm3, nEl=1);
    G4Element *elPb = new G4Element(name="Lead", symbol="Pb", iz=82., a=207.19*g/mole);
    thePb->AddElement( elPb, 1 );
    G4Material *theZr = new G4Material(name="Zirconium", density=1.53*g/cm3, nEl=1); 
    G4Element *elZr = new G4Element(name="Zirconium", symbol="Zr", iz=40., a=91.224*g/mole);
    theZr->AddElement( elZr, 1 );
    
    G4int numberOfMaterials = 2;
    G4Material *theMaterials[2];
    theMaterials[0] = thePb;
    theMaterials[1] = theZr;

    // ----------- here all materials have been defined -----------
    
    //G4Element::DumpInfo(); 
    //G4Material::DumpInfo();
    
    // ----------- the following is needed for building a track ------------
    
    static const G4MaterialTable *theMaterialTable = G4Material::GetMaterialTable();
    G4int imat = 0;
    G4Box* theFrame = new G4Box ( "Frame",10*m, 10*m, 10*m );
    
    G4LogicalVolume *LogicalFrame = new G4LogicalVolume( theFrame,
                                                         (*theMaterialTable)(imat),
                                                         "LFrame", 0, 0, 0 );
    
    G4PVPlacement *PhysicalFrame = new G4PVPlacement( 0, G4ThreeVector(),
                                                     "PFrame", LogicalFrame, 0, false, 0 );
    
    G4RotationMatrix theNull;
    G4ThreeVector theCenter(0,0,0);
    G4GRSVolume * theTouchable = new G4GRSVolume(PhysicalFrame, &theNull, theCenter);
    // ----------- now get all particles of interest ---------
   
    G4BosonConstructor Bosons;
    Bosons.ConstructParticle();

    G4MesonConstructor Mesons;
    Mesons.ConstructParticle();

    G4LeptonConstructor Leptons;
    Leptons.ConstructParticle();

    G4BaryonConstructor Baryons;
    Baryons.ConstructParticle();
 
    G4int numberOfParticles = 25;
    G4ParticleDefinition *theParticles[25];
    G4ParticleDefinition *theProton = G4Proton::ProtonDefinition();
    theProton->SetCuts(0.01*mm);
    theParticles[ 0] = theProton;
    G4ParticleDefinition *theAntiProton = G4AntiProton::AntiProtonDefinition();
    theParticles[ 1] = theAntiProton;
    G4ParticleDefinition *theNeutron = G4Neutron::NeutronDefinition();
    theParticles[ 2] = theNeutron;
    G4ParticleDefinition *theAntiNeutron = G4AntiNeutron::AntiNeutronDefinition();
    theParticles[ 3] = theAntiNeutron;
    G4ParticleDefinition *thePionPlus = G4PionPlus::PionPlusDefinition();
    theParticles[ 4] = thePionPlus;
    G4ParticleDefinition *thePionMinus = G4PionMinus::PionMinusDefinition();
    theParticles[ 5] = thePionMinus;
    G4ParticleDefinition *theDeuteron = G4Deuteron::DeuteronDefinition();
    theParticles[ 6] = theDeuteron;
    G4ParticleDefinition *theTriton = G4Triton::TritonDefinition();
    theParticles[ 7] = theTriton;
    G4ParticleDefinition *theAlpha = G4Alpha::AlphaDefinition();
    theParticles[ 8] = theAlpha;
    G4ParticleDefinition *theKaonPlus = G4KaonPlus::KaonPlusDefinition();
    theParticles[ 9] = theKaonPlus;
    G4ParticleDefinition *theKaonMinus = G4KaonMinus::KaonMinusDefinition();
    theParticles[10] = theKaonMinus;
    G4ParticleDefinition *theKaonLong = G4KaonZeroLong::KaonZeroLongDefinition();
    theParticles[11] = theKaonLong;
    G4ParticleDefinition *theKaonShort = G4KaonZeroShort::KaonZeroShortDefinition();
    theParticles[12] = theKaonShort;
    G4ParticleDefinition *theLambda = G4Lambda::LambdaDefinition();
    theParticles[13] = theLambda;
    G4ParticleDefinition *theAntiLambda = G4AntiLambda::AntiLambdaDefinition();
    theParticles[14] = theAntiLambda;
    G4ParticleDefinition *theSigmaPlus = G4SigmaPlus::SigmaPlusDefinition();
    theParticles[15] = theSigmaPlus;
    G4ParticleDefinition *theAntiSigmaPlus = G4AntiSigmaPlus::AntiSigmaPlusDefinition();
    theParticles[16] = theAntiSigmaPlus;
    G4ParticleDefinition *theSigmaMinus = G4SigmaMinus::SigmaMinusDefinition();
    theParticles[17] = theSigmaMinus;
    G4ParticleDefinition *theAntiSigmaMinus = G4AntiSigmaMinus::AntiSigmaMinusDefinition();
    theParticles[18] = theAntiSigmaMinus;
    G4ParticleDefinition *theXiMinus = G4XiMinus::XiMinusDefinition();
    theParticles[19] = theXiMinus;
    G4ParticleDefinition *theAntiXiMinus = G4AntiXiMinus::AntiXiMinusDefinition();
    theParticles[20] = theAntiXiMinus;
    G4ParticleDefinition *theXiZero = G4XiZero::XiZeroDefinition();
    theParticles[21] = theXiZero;
    G4ParticleDefinition *theAntiXiZero = G4AntiXiZero::AntiXiZeroDefinition();
    theParticles[22] = theAntiXiZero;
    G4ParticleDefinition *theOmegaMinus = G4OmegaMinus::OmegaMinusDefinition();
    theParticles[23] = theOmegaMinus;
    G4ParticleDefinition *theAntiOmegaMinus = G4AntiOmegaMinus::AntiOmegaMinusDefinition();
    theParticles[24] = theAntiOmegaMinus;
    
    //------ all the particles are Done ----------
    //------ Processes definitions Follow ---------

    // this will be the model class for low energies
//    G4HardronNuclearChips * theTheoModel = new G4HardronNuclearChips;
    G4Evaporation * theEvaporation = new G4Evaporation;
    G4FermiBreakUp * theFermiBreakUp = new G4FermiBreakUp;
    G4StatMF * theMF = new G4StatMF;
    G4ExcitationHandler * theHandler = new G4ExcitationHandler;
        theHandler->SetEvaporation(theEvaporation);
        theHandler->SetFermiModel(theFermiBreakUp);
        theHandler->SetMultiFragmentation(theMF);
        theHandler->SetMaxAandZForFermiBreakUp(12, 6);
        theHandler->SetMinEForMultiFrag(3*MeV);
	
    G4PreCompoundModel * theTheoModel = new G4PreCompoundModel(theHandler);
       
    
    G4HadronInelasticProcess *theProcesses[25];
    G4ProcessManager *theProtonProcessManager = new G4ProcessManager(theProton);
    theProton->SetProcessManager(theProtonProcessManager);

      //      theProton->GetProcessManager();
    G4ProtonInelasticProcess *theProtonInelasticProcess =
      new G4ProtonInelasticProcess();
//    theProtonInelasticProcess->RegisterMe( theParaModel );
    theProtonInelasticProcess->RegisterMe( theTheoModel );
    theProtonProcessManager->AddDiscreteProcess( theProtonInelasticProcess );
    theProcesses[0] = theProtonInelasticProcess;
   
    //    G4ProcessManager *theAntiProtonProcessManager =
    //      theAntiProton->GetProcessManager();

    G4ProcessManager *theAntiProtonProcessManager = new G4ProcessManager(theAntiProton);
    theAntiProton->SetProcessManager(theAntiProtonProcessManager);

    G4AntiProtonInelasticProcess *theAntiProtonInelasticProcess =
      new G4AntiProtonInelasticProcess(); 
    theAntiProtonInelasticProcess->RegisterMe( theTheoModel );
    theAntiProtonProcessManager->AddDiscreteProcess( theAntiProtonInelasticProcess );
    theProcesses[1] = theAntiProtonInelasticProcess;
    
    //    G4ProcessManager *theNeutronProcessManager =
    //      theNeutron->GetProcessManager();
    G4ProcessManager *theNeutronProcessManager = new G4ProcessManager(theNeutron);
    theNeutron->SetProcessManager(theNeutronProcessManager);
    G4NeutronInelasticProcess *theNeutronInelasticProcess =
      new G4NeutronInelasticProcess(); 
    theNeutronInelasticProcess->RegisterMe( theTheoModel );
    theNeutronProcessManager->AddDiscreteProcess( theNeutronInelasticProcess );
    theProcesses[2] = theNeutronInelasticProcess;
    
    //    G4ProcessManager *theAntiNeutronProcessManager =
    //      theAntiNeutron->GetProcessManager();

    G4ProcessManager *theAntiNeutronProcessManager = new G4ProcessManager(theAntiNeutron);
    theAntiNeutron->SetProcessManager(theAntiNeutronProcessManager);
    G4AntiNeutronInelasticProcess *theAntiNeutronInelasticProcess =
      new G4AntiNeutronInelasticProcess();
    theAntiNeutronInelasticProcess->RegisterMe( theTheoModel );
    theAntiNeutronProcessManager->AddDiscreteProcess( theAntiNeutronInelasticProcess );
    theProcesses[3] = theAntiNeutronInelasticProcess;

    //    G4ProcessManager *thePionPlusProcessManager =
    //      thePionPlus->GetProcessManager();
    G4ProcessManager *thePionPlusProcessManager = new G4ProcessManager(thePionPlus);
    thePionPlus->SetProcessManager(thePionPlusProcessManager);
    G4PionPlusInelasticProcess *thePionPlusInelasticProcess =
      new G4PionPlusInelasticProcess();
    thePionPlusInelasticProcess->RegisterMe( theTheoModel );
    thePionPlusProcessManager->AddDiscreteProcess( thePionPlusInelasticProcess );
    theProcesses[4] = thePionPlusInelasticProcess;
    
    //    G4ProcessManager *thePionMinusProcessManager =
    //      thePionMinus->GetProcessManager();
    G4ProcessManager *thePionMinusProcessManager = new G4ProcessManager(thePionMinus);
    thePionMinus->SetProcessManager(thePionMinusProcessManager);
    G4PionMinusInelasticProcess *thePionMinusInelasticProcess =
      new G4PionMinusInelasticProcess();
    thePionMinusInelasticProcess->RegisterMe( theTheoModel );
    thePionMinusProcessManager->AddDiscreteProcess( thePionMinusInelasticProcess );
    theProcesses[5] = thePionMinusInelasticProcess;
    
    //    G4ProcessManager *theDeuteronProcessManager =
    //      theDeuteron->GetProcessManager();
    G4ProcessManager *theDeuteronProcessManager = new G4ProcessManager(theDeuteron);
    theDeuteron->SetProcessManager(theDeuteronProcessManager);
    G4DeuteronInelasticProcess *theDeuteronInelasticProcess =
      new G4DeuteronInelasticProcess();
    theDeuteronInelasticProcess->RegisterMe( theTheoModel );
    theDeuteronProcessManager->AddDiscreteProcess( theDeuteronInelasticProcess );
    theProcesses[6] = theDeuteronInelasticProcess;
    
    //    G4ProcessManager *theTritonProcessManager =
    //      theTriton->GetProcessManager();
    G4ProcessManager *theTritonProcessManager = new G4ProcessManager(theTriton);
    theTriton->SetProcessManager(theTritonProcessManager);
    G4TritonInelasticProcess *theTritonInelasticProcess =
      new G4TritonInelasticProcess();
    theTritonInelasticProcess->RegisterMe( theTheoModel );
    theTritonProcessManager->AddDiscreteProcess( theTritonInelasticProcess );
    theProcesses[7] = theTritonInelasticProcess;
    
    //    G4ProcessManager *theAlphaProcessManager =
    //      theAlpha->GetProcessManager();
    G4ProcessManager *theAlphaProcessManager = new G4ProcessManager(theAlpha);
    theAlpha->SetProcessManager(theAlphaProcessManager);
    G4AlphaInelasticProcess *theAlphaInelasticProcess =
      new G4AlphaInelasticProcess();
    theAlphaInelasticProcess->RegisterMe( theTheoModel );
    theAlphaProcessManager->AddDiscreteProcess( theAlphaInelasticProcess );
    theProcesses[8] = theAlphaInelasticProcess;
    
    //    G4ProcessManager *theKaonPlusProcessManager =
    //      theKaonPlus->GetProcessManager();
    G4ProcessManager *theKaonPlusProcessManager = new G4ProcessManager(theKaonPlus);
    theKaonPlus->SetProcessManager(theKaonPlusProcessManager);
    G4KaonPlusInelasticProcess *theKaonPlusInelasticProcess =
      new G4KaonPlusInelasticProcess();
    theKaonPlusInelasticProcess->RegisterMe( theTheoModel );
    theKaonPlusProcessManager->AddDiscreteProcess( theKaonPlusInelasticProcess );
    theProcesses[9] = theKaonPlusInelasticProcess;
    
    //    G4ProcessManager *theKaonMinusProcessManager =
    //      theKaonMinus->GetProcessManager();
    G4ProcessManager *theKaonMinusProcessManager = new G4ProcessManager(theKaonMinus);
    theKaonMinus->SetProcessManager(theKaonMinusProcessManager);
    G4KaonMinusInelasticProcess *theKaonMinusInelasticProcess =
      new G4KaonMinusInelasticProcess();
    theKaonMinusInelasticProcess->RegisterMe( theTheoModel );
    theKaonMinusProcessManager->AddDiscreteProcess( theKaonMinusInelasticProcess );
    theProcesses[10] = theKaonMinusInelasticProcess;

    //    G4ProcessManager *theKaonLongProcessManager =
    //      theKaonLong->GetProcessManager();
    G4ProcessManager *theKaonLongProcessManager = new G4ProcessManager(theKaonLong);
    theKaonLong->SetProcessManager(theKaonLongProcessManager);
    G4KaonZeroLInelasticProcess *theKaonLongInelasticProcess =
      new G4KaonZeroLInelasticProcess();
    theKaonLongInelasticProcess->RegisterMe( theTheoModel );
    theKaonLongProcessManager->AddDiscreteProcess( theKaonLongInelasticProcess );
    theProcesses[11] = theKaonLongInelasticProcess;
    
    //    G4ProcessManager *theKaonShortProcessManager =
    //      theKaonShort->GetProcessManager();
    G4ProcessManager *theKaonShortProcessManager = new G4ProcessManager(theKaonShort);
    theKaonShort->SetProcessManager(theKaonShortProcessManager);
    G4KaonZeroSInelasticProcess *theKaonShortInelasticProcess =
      new G4KaonZeroSInelasticProcess();
    theKaonShortInelasticProcess->RegisterMe( theTheoModel );
    theKaonShortProcessManager->AddDiscreteProcess( theKaonShortInelasticProcess );
    theProcesses[12] = theKaonShortInelasticProcess;
    
    //    G4ProcessManager *theLambdaProcessManager =
    //      theLambda->GetProcessManager();
    G4ProcessManager *theLambdaProcessManager = new G4ProcessManager(theLambda);
    theLambda->SetProcessManager(theLambdaProcessManager);
    G4LambdaInelasticProcess *theLambdaInelasticProcess =
      new G4LambdaInelasticProcess();
    theLambdaInelasticProcess->RegisterMe( theTheoModel );
    theLambdaProcessManager->AddDiscreteProcess( theLambdaInelasticProcess );
    theProcesses[13] = theLambdaInelasticProcess;
    
    //    G4ProcessManager *theAntiLambdaProcessManager =
    //      theAntiLambda->GetProcessManager();
    G4ProcessManager *theAntiLambdaProcessManager = new G4ProcessManager(theAntiLambda);
    theAntiLambda->SetProcessManager(theAntiLambdaProcessManager);
    G4AntiLambdaInelasticProcess *theAntiLambdaInelasticProcess =
      new G4AntiLambdaInelasticProcess();
    theAntiLambdaInelasticProcess->RegisterMe( theTheoModel );
    theAntiLambdaProcessManager->AddDiscreteProcess( theAntiLambdaInelasticProcess );
    theProcesses[14] = theAntiLambdaInelasticProcess;
    
    //    G4ProcessManager *theSigmaPlusProcessManager =
    //      theSigmaPlus->GetProcessManager();
    G4ProcessManager *theSigmaPlusProcessManager = new G4ProcessManager(theSigmaPlus);
    theSigmaPlus->SetProcessManager(theSigmaPlusProcessManager);
    G4SigmaPlusInelasticProcess *theSigmaPlusInelasticProcess =
       new G4SigmaPlusInelasticProcess();
    theSigmaPlusInelasticProcess->RegisterMe( theTheoModel );
    theSigmaPlusProcessManager->AddDiscreteProcess( theSigmaPlusInelasticProcess );
    theProcesses[15] = theSigmaPlusInelasticProcess;
    
    //    G4ProcessManager *theAntiSigmaPlusProcessManager =
    //      theAntiSigmaPlus->GetProcessManager();
    G4ProcessManager *theAntiSigmaPlusProcessManager = new G4ProcessManager(theAntiSigmaPlus);
    theAntiSigmaPlus->SetProcessManager(theAntiSigmaPlusProcessManager);
    G4AntiSigmaPlusInelasticProcess *theAntiSigmaPlusInelasticProcess =
      new G4AntiSigmaPlusInelasticProcess();
    theAntiSigmaPlusInelasticProcess->RegisterMe( theTheoModel );
    theAntiSigmaPlusProcessManager->AddDiscreteProcess( theAntiSigmaPlusInelasticProcess );
    theProcesses[16] = theAntiSigmaPlusInelasticProcess;
    
    //    G4ProcessManager *theSigmaMinusProcessManager =
    //      theSigmaMinus->GetProcessManager();
    G4ProcessManager *theSigmaMinusProcessManager = new G4ProcessManager(theSigmaMinus);
    theSigmaMinus->SetProcessManager(theSigmaMinusProcessManager);
    G4SigmaMinusInelasticProcess *theSigmaMinusInelasticProcess =
      new G4SigmaMinusInelasticProcess();
    theSigmaMinusInelasticProcess->RegisterMe( theTheoModel );
    theSigmaMinusProcessManager->AddDiscreteProcess( theSigmaMinusInelasticProcess );
    theProcesses[17] = theSigmaMinusInelasticProcess;
    
    //    G4ProcessManager *theAntiSigmaMinusProcessManager =
    //      theAntiSigmaMinus->GetProcessManager();
    G4ProcessManager *theAntiSigmaMinusProcessManager = new G4ProcessManager(theAntiSigmaMinus);
    theAntiSigmaMinus->SetProcessManager(theAntiSigmaMinusProcessManager);
    G4AntiSigmaMinusInelasticProcess *theAntiSigmaMinusInelasticProcess =
      new G4AntiSigmaMinusInelasticProcess();
    theAntiSigmaMinusInelasticProcess->RegisterMe( theTheoModel );
    theAntiSigmaMinusProcessManager->AddDiscreteProcess( theAntiSigmaMinusInelasticProcess );
    theProcesses[18] = theAntiSigmaMinusInelasticProcess;
    
    //    G4ProcessManager *theXiMinusProcessManager =
    //      theXiMinus->GetProcessManager();
    G4ProcessManager *theXiMinusProcessManager = new G4ProcessManager(theXiMinus);
    theXiMinus->SetProcessManager(theXiMinusProcessManager);
    G4XiMinusInelasticProcess *theXiMinusInelasticProcess =
      new G4XiMinusInelasticProcess();
    theXiMinusInelasticProcess->RegisterMe( theTheoModel );
    theXiMinusProcessManager->AddDiscreteProcess( theXiMinusInelasticProcess );
    theProcesses[19] = theXiMinusInelasticProcess;
    
    //    G4ProcessManager *theAntiXiMinusProcessManager =
    //      theAntiXiMinus->GetProcessManager();
    G4ProcessManager *theAntiXiMinusProcessManager = new G4ProcessManager(theAntiXiMinus);
    theAntiXiMinus->SetProcessManager(theAntiXiMinusProcessManager);
    G4AntiXiMinusInelasticProcess *theAntiXiMinusInelasticProcess =
      new G4AntiXiMinusInelasticProcess();
    theAntiXiMinusInelasticProcess->RegisterMe( theTheoModel );
    theAntiXiMinusProcessManager->AddDiscreteProcess( theAntiXiMinusInelasticProcess );
    theProcesses[20] = theAntiXiMinusInelasticProcess;
    
    //    G4ProcessManager *theXiZeroProcessManager =
    //      theXiZero->GetProcessManager();
    G4ProcessManager *theXiZeroProcessManager = new G4ProcessManager(theXiZero);
    theXiZero->SetProcessManager(theXiZeroProcessManager);
    G4XiZeroInelasticProcess *theXiZeroInelasticProcess =
      new G4XiZeroInelasticProcess();
    theXiZeroInelasticProcess->RegisterMe( theTheoModel );
    theXiZeroProcessManager->AddDiscreteProcess( theXiZeroInelasticProcess );
    theProcesses[21] = theXiZeroInelasticProcess;
    
    //    G4ProcessManager *theAntiXiZeroProcessManager =
    //      theAntiXiZero->GetProcessManager();
    G4ProcessManager *theAntiXiZeroProcessManager = new G4ProcessManager(theAntiXiZero);
    theAntiXiZero->SetProcessManager(theAntiXiZeroProcessManager);
    G4AntiXiZeroInelasticProcess *theAntiXiZeroInelasticProcess =
      new G4AntiXiZeroInelasticProcess();
    theAntiXiZeroInelasticProcess->RegisterMe( theTheoModel );
    theAntiXiZeroProcessManager->AddDiscreteProcess( theAntiXiZeroInelasticProcess );
    theProcesses[22] = theAntiXiZeroInelasticProcess;
    
    //    G4ProcessManager *theOmegaMinusProcessManager =
    //      theOmegaMinus->GetProcessManager();
    G4ProcessManager *theOmegaMinusProcessManager = new G4ProcessManager(theOmegaMinus);
    theOmegaMinus->SetProcessManager(theOmegaMinusProcessManager);
    G4OmegaMinusInelasticProcess *theOmegaMinusInelasticProcess =
      new G4OmegaMinusInelasticProcess();
    theOmegaMinusInelasticProcess->RegisterMe( theTheoModel );
    theOmegaMinusProcessManager->AddDiscreteProcess( theOmegaMinusInelasticProcess );
    theProcesses[23] = theOmegaMinusInelasticProcess;
    
    //    G4ProcessManager *theAntiOmegaMinusProcessManager =
    //      theAntiOmegaMinus->GetProcessManager();
    G4ProcessManager *theAntiOmegaMinusProcessManager = new G4ProcessManager(theAntiOmegaMinus);
    theAntiOmegaMinus->SetProcessManager(theAntiOmegaMinusProcessManager);
    G4AntiOmegaMinusInelasticProcess *theAntiOmegaMinusInelasticProcess =
      new G4AntiOmegaMinusInelasticProcess();
    theAntiOmegaMinusInelasticProcess->RegisterMe( theTheoModel );
    theAntiOmegaMinusProcessManager->AddDiscreteProcess( theAntiOmegaMinusInelasticProcess );
    theProcesses[24] = theAntiOmegaMinusInelasticProcess;
    
    G4ForceCondition *condition = new G4ForceCondition;
    *condition = NotForced;
    
    G4ParticleMomentum theDirection( 0., 0., 1. );
    G4ThreeVector aPosition( 0., 0., 0. );
    G4double aTime = 0.0;
    
    G4StepPoint aStepPoint;
    G4Step aStep;
    aStep.SetPreStepPoint(&aStepPoint);
    G4double meanFreePath;
    G4double incomingEnergy;
    G4ParticleDefinition *currentDefinition;
    G4Track *currentTrack;
    
    G4int kl, kr;
    do {
      G4cout << " 0) Lead      1)  Zirconium"<< G4endl;
      G4cout << "Please enter the material code" << G4endl;
      G4cout << "\tFrom: " << G4std::flush;
      G4cin >> kl;
      G4cout << "\tTo: " << G4std::flush;
      G4cin >> kr;
    } while( kl < 0 || kr >= numberOfMaterials || kl > kr );
    
    G4int il, ir;
    do {
      G4cout << " 0) proton        1) anti_proton    2) neutron       3) anti_neutron" << G4endl;
      G4cout << " 4) pi+           5) pi-            6) deuteron      7) triton" << G4endl;
      G4cout << " 8) alpha         9) kaon+         10) kaon-        11) kaon0L" << G4endl;
      G4cout << "12) kaon0S       13) lambda        14) anti_lambda  15) sigma+" << G4endl;
      G4cout << "16) anti_sigma+  17) sigma-        18) anti_sigma-  19) xi-" << G4endl;
      G4cout << "20) anti_xi-     21) xi0           22) anti_xi0     23) omega-" << G4endl;
      G4cout << "24) anti_omega-" << G4endl;
      G4cout << "Please enter the particle type code"<< G4endl;
      G4cout << "\tFrom: " << G4std::flush;
      G4cin >> il;
      G4cout << "\tTo: " << G4std::flush;
      G4cin >> ir;
    } while( il < 0 || ir > 24 || il > ir );

    G4cout << "Please enter the initial kinetic energy (GeV): "<< G4std::flush;
    G4cin >> incomingEnergy;
    incomingEnergy *= GeV;

    
    G4int nEvents;
    G4cout << "Please enter the number of events: "<< G4std::flush;
    G4cin >> nEvents;
    
    if( kl == kr && il == ir ) {
    }
   
    G4int nCPU = 0;
    G4double dCPU = 0.0;
    G4Timer timerEvent;
    G4Timer timerTotal;
    timerTotal.Start();
    
    // Prepare the analysis
    G4String file("../logs/liste.");
    G4String fileName; // p_pb
    G4String eString; // .80MeV
    G4double weight; // sigma
    G4ProtonInelasticCrossSection theCrossSections;
    maximum = incomingEnergy; 
    G4DynamicParticle * aDummyParticle = new G4DynamicParticle(G4Proton::ProtonDefinition(), G4ParticleMomentum(1.,0.,0.), incomingEnergy);
    if(kl == 0) 
    {
      fileName = "p_pb";
      weight = theCrossSections.GetCrossSection(aDummyParticle, elPb, 273*kelvin);
    }
    if(kl == 1) 
    {
      fileName = "p_zr";
      weight = theCrossSections.GetCrossSection(aDummyParticle, elZr, 273*kelvin);
    }
    if(abs(incomingEnergy-113*MeV)<1*MeV) eString = ".113mev";
    if(abs(incomingEnergy-120*MeV)<1*MeV) eString = ".120mev";
    if(abs(incomingEnergy-160*MeV)<1*MeV) eString = ".160mev";
    if(abs(incomingEnergy-256*MeV)<1*MeV) eString = ".256mev";
    if(abs(incomingEnergy-25*MeV)<1*MeV) eString = ".25mev";
    if(abs(incomingEnergy-45*MeV)<1*MeV) eString = ".45mev";
    if(abs(incomingEnergy-597*MeV)<1*MeV) eString = ".597mev";
    if(abs(incomingEnergy-800*MeV)<1*MeV) eString = ".800mev";
    if(abs(incomingEnergy-80*MeV)<1*MeV) eString = ".80mev";
    fileName = fileName+eString;
    G4cout << "running on"<<fileName<<G4endl;
    G4String it;
    it = file+fileName;
    ANAParticleInfo theInformation(weight, it);
 
    G4int k, i;
    for( k = kl; k <= kr; k++ )
      for( i = il; i <= ir; i++ ) {
	theParticles[i]->SetCuts( (i ? 1.0 : 0.01)*mm );
	LogicalFrame->SetMaterial( theMaterials[k] ); 
    
    // --------- Test the PostStepDoIt now  --------------
    
	G4ParticleChange *aFinalState;
	LogicalFrame->SetMaterial( theMaterials[k] ); 
	G4DynamicParticle aParticle( theParticles[i], theDirection, incomingEnergy );
	G4double currentPx, currentPy, currentPz;
	G4double kineticEnergy, totalEnergy, currentMass;
    
	G4int eventCounter = 0;
	G4int l;
	for( l=0; l<nEvents; ++l ) {
	  aParticle.SetDefinition( theParticles[i] );
	  aParticle.SetMomentumDirection( theDirection );
	  aParticle.SetKineticEnergy( incomingEnergy );
      
	  G4Track *aTrack = new G4Track( &aParticle, aTime, aPosition );
	  aTrack->SetTouchable( theTouchable );
	  aTrack->SetStep( &aStep );
	  aStep.SetTrack( aTrack );
          aStepPoint.SetTouchable(theTouchable);
	  aStepPoint.SetMaterial(theMaterials[k]);
          aStep.SetPreStepPoint(&aStepPoint);
	  aStep.SetPostStepPoint(&aStepPoint);
	  G4cout << "Event:" << G4std::setw(4) << l+1 << 
	            " Material:" << theMaterials[k]->GetName() <<
	            " Particle:" << theParticles[i]->GetParticleName() << '\t' << G4endl;

	  G4cerr << "Next event is: #"<<l<<G4endl;
	  timerEvent.Start();
	  aFinalState = (G4ParticleChange * ) (theProcesses[i]->PostStepDoIt( *aTrack, aStep ) );
	  timerEvent.Stop();

	  G4int nSec = aFinalState->GetNumberOfSecondaries();
      
	  // prepare the analysis current quantities
      
          G4cout << "NUMBER OF SECONDARIES="<<aFinalState->GetNumberOfSecondaries()<<G4endl;
	  G4int istart, ipar;  
	  double etWithinCuts = 0;    
	  G4int ii;  
          G4Track * second;
	  G4DynamicParticle * aSec;
          G4int isec;
          for(isec=0;isec<aFinalState->GetNumberOfSecondaries();isec++)
          {
            second = aFinalState->GetSecondary(isec);
            aSec = second->GetDynamicParticle();
	    G4cout << "SECONDARIES info";
            G4cout << aSec->GetDefinition()->GetPDGCharge()<<" ";
	    if(aSec->GetDefinition()->GetPDGEncoding())
              G4cout << aSec->GetDefinition()->GetPDGEncoding()<<" ";
	    else
	      G4cout << aSec->GetDefinition()->GetPDGCharge()*1000+aSec->GetDefinition()->GetBaryonNumber()
	      <<" ";
            G4cout << aSec->GetTotalEnergy()<<" ";
            G4cout << aSec->GetMomentum();
	    G4cout << (1-isec)*aFinalState->GetNumberOfSecondaries();
	    G4cout << G4endl;
	    // analysis part;
	    G4int theCHIPSCode = aSec->GetDefinition()->GetPDGEncoding();
	    if(theCHIPSCode == 0)
	    {
	      theCHIPSCode =
	      aSec->GetDefinition()->GetPDGCharge()*1000+aSec->GetDefinition()->GetBaryonNumber();
	    }
	    ANAParticle aPart(theCHIPSCode,
	                      aSec->GetMomentum().x(),
			      aSec->GetMomentum().y(),
			      aSec->GetMomentum().z(),
			      aSec->GetTotalEnergy());
            theInformation.ProcessOne(aPart);
	    delete aSec;
	    delete second;
          }
          delete aTrack;
	  aFinalState->Clear();
	  if(l==100*(l/100)) 
	  {
	    theInformation.Analyse();
	    theInformation.Plot(fileName, l);
	  }
	} // end of event loop
      } // end of i loop
    timerTotal.Stop();    
    G4cout << "Terminating successfully"<<endl;
    G4cout << "\tTotal time " << timerTotal <<
              "\n\tPure function time:" << dCPU <<
              "\n\tTotal number of events:" << nCPU << G4endl; 
    return EXIT_SUCCESS;
}
 /* end of code */
 
