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
// $Id: InelasticAlphaTest.cc,v 1.2 2003-10-08 12:19:48 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Johannes Peter Wellisch, 22.Apr 1997: full test-suite coded.    
#include "G4ios.hh"
#include <fstream>
#include <iomanip>
 
#include "G4Material.hh"
 
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
 
#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4StatEvaporation.hh"
#include "G4StatCompetitiveFission.hh"
#include "G4StatFermiBreakUp.hh"
#include "G4DummyMF.hh"
#include "G4HadronKineticModel.hh"
#include "G4Annihilator.hh"
#include "G4ElasticScatterer.hh"
#include "G4InElasticScatterer.hh"
#include "G4CrossSectionBase.hh"
#include "G4RKFieldIntegrator.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4CrossSectionBase.hh"
#include "G4LEProtonInelastic.hh"
#include "G4StringModel.hh"
#include "G4PreCompoundModel.hh"
#include "G4RKFieldIntegrator.hh"
#include "G4FTFModel.hh"

#include "G4DynamicParticle.hh"
#include "G4Proton.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4GRSVolume.hh"
#include "G4Step.hh"

#include "g4templates.hh"
 
#include "NametoGeant3Number.cc"
 
 int main()
  {
    G4cout.setf( std::ios::scientific, std::ios::floatfield );
    std::ofstream outFile( "InelasticAlpha.listing.GetMeanFreePath", std::ios::out);
    outFile.setf( std::ios::scientific, std::ios::floatfield );
    std::ofstream outFile1( "InelasticAlpha.listing.DoIt", std::ios::out);
    outFile1.setf( std::ios::scientific, std::ios::floatfield );
    
    FILE *g4data;
    if( (g4data = fopen("g4data","w")) == NULL )
    {
      G4cerr << "cannot open g4data" << G4endl;
      exit( EXIT_FAILURE );
    }
    RanecuEngine theEngine;
    HepRandom::setTheEngine( &theEngine );
    theEngine.showStatus();
    G4cout << G4endl;
    
    G4String name, symbol;
    G4double a, iz, z, density;
    G4int nEl;
    
    G4Material *theCu = new G4Material(name="Copper", density=8.96*g/cm3, nEl=1);
    G4Element *elCu = new G4Element(name="Copper", symbol="Cu", iz=29., a=63.55*g/mole);
    theCu->AddElement( elCu, 1 );
    G4Material *thePb = new G4Material(name="Lead", density=11.35*g/cm3, nEl=1);
    G4Element *elPb = new G4Element(name="Lead", symbol="Pb", iz=82., a=207.19*g/mole);
    thePb->AddElement( elPb, 1 );
    G4Material *theFe = new G4Material(name="Iron", density=7.87*g/cm3, nEl=1);
    G4Element *elFe = new G4Element(name="Iron", symbol="Fe", iz=26., a=55.85*g/mole);
    theFe->AddElement( elFe, 1 );
    G4Material *theW  = new G4Material(name="Tungsten", density=19.30*g/cm3, nEl=1);
    G4Element *elW  = new G4Element(name="Tungston", symbol="W", iz=74., a=183.85*g/mole);
    theW->AddElement( elW, 1 );
    G4Material *theLAr= new G4Material(name="LArgon", density=1.393*g/cm3, nEl=1);
    G4Element *elAr  = new G4Element(name="Argon", symbol="Ar", iz=18., a=39.95*g/mole);
    theLAr->AddElement( elAr, 1 );
    G4Material *thePS = new G4Material(name="PolyStyrene", density=1.032*g/cm3, nEl=2);
    G4Element *elC = new G4Element(name="Carbon", symbol="C", iz=6., a=12.01*g/mole);
    G4Element *elH = new G4Element(name="Hydrogen", symbol="H", iz=1., a=1.01*g/mole);
    thePS->AddElement( elC, 8 );
    thePS->AddElement( elH, 8 );
    G4Material *thePbWO4 = new G4Material(name="PbWO4", density=12.0*g/cm3, nEl=3);
    // approximative number
    G4Element *elO = new G4Element(name="Oxygen", symbol="O", iz=8., a=15.9994*g/mole);
    thePbWO4->AddElement( elPb, 1 );
    thePbWO4->AddElement( elW,  1 );
    thePbWO4->AddElement( elO,  4 );
    // approximate numbers for O
    G4Material *theO = new G4Material(name="Oxygen", density=1.1*g/cm3, nEl=1);
    theO->AddElement( elO,  1 );
    G4Material *theBe = new G4Material(name="Beryllium", density=1.848*g/cm3, nEl=1);
    G4Element *elBe  = new G4Element(name="Beryllium", symbol="Be", iz=4., a=9.01*g/mole);
    theBe->AddElement( elBe, 1 );
    G4Material *theAl = new G4Material(name="Aluminium", density=2.70*g/cm3, nEl=1);
    G4Element *elAl  = new G4Element(name="Aluminium", symbol="Al", iz=13., a=26.98*g/mole);
    theAl->AddElement( elAl, 1 );
    G4Material *theU = new G4Material(name="Uranium", density=18.95*g/cm3, nEl=1);
    G4Element *elU  = new G4Element(name="Uranium", symbol="U", iz=92., a=238.03*g/mole);
    theU->AddElement( elU, 1 );
    G4Material *theBGO = new G4Material(name="BGO", density=2.15*g/cm3, nEl=3);
    G4Element *elBi = new G4Element(name="Bismuth", symbol="Bi", iz=83., a=208.98*g/mole);
    G4Element *elGe = new G4Element(name="Germanium", symbol="Ge", iz=32., a=72.59*g/mole);
    theBGO->AddElement( elBi, 4 );
    theBGO->AddElement( elGe, 3 );
    theBGO->AddElement( elO, 12 );
    G4Material *theNaI = new G4Material(name="NaI", density=3.67*g/cm3, nEl=2);
    G4Element *elNa = new G4Element(name="Sodium", symbol="Na", iz=11., a=22.990*g/mole);
    G4Element *elI = new G4Element(name="Iodine", symbol="I", iz=53., a=126.904*g/mole);
    theNaI->AddElement( elNa, 1 );
    theNaI->AddElement( elI, 1 );
    G4Material *theCsI = new G4Material(name="CsI", density=4.53*g/cm3, nEl=2);
    G4Element *elCs = new G4Element(name="Cesium", symbol="Cs", iz=55., a=132.905*g/mole);
    theCsI->AddElement( elCs, 1 );
    theCsI->AddElement( elI, 1 );
    G4Material *theKapton = new G4Material(name="Kapton", density=1.53*g/cm3, nEl=4); 
    // formula: private communications, see mail.
    theKapton->AddElement( elC, 22 );
    theKapton->AddElement( elH, 10 );
    theKapton->AddElement( elO,  5 );
    G4Element *elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a=14.007*g/mole);
    theKapton->AddElement( elN, 2 );
    
    G4int numberOfMaterials = 15;
    G4Material *theMaterials[15];
    theMaterials[ 0] = theCu;
    theMaterials[ 1] = thePb;
    theMaterials[ 2] = theFe;
    theMaterials[ 3] = theW;
    theMaterials[ 4] = theLAr;
    theMaterials[ 5] = thePS;
    theMaterials[ 6] = thePbWO4;
    theMaterials[ 7] = theO;
    theMaterials[ 8] = theBe;
    theMaterials[ 9] = theAl;
    theMaterials[10] = theU;
    theMaterials[11] = theBGO;
    theMaterials[12] = theNaI;
    theMaterials[13] = theCsI;
    theMaterials[14] = theKapton;
    
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
    
    G4int numberOfParticles = 25;
    G4ParticleDefinition *theParticles[25];
    G4ParticleDefinition *theProton = G4Proton::ProtonDefinition();
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

    // this will be the model class for high energies
    G4TheoFSGenerator * theTheoModel = new G4TheoFSGenerator;
    
    // all models for treatment of thermal nucleus 
    G4StatCompetitiveFission * theCompFission = new G4StatCompetitiveFission;
    G4StatEvaporation * theEvaporation = new G4StatEvaporation;
        theEvaporation->SetFission(theCompFission);
    G4StatFermiBreakUp * theFermiBreakUp = new G4StatFermiBreakUp;
    G4DummyMF * theMF = new G4DummyMF;

    // Evaporation logic
    G4ExcitationHandler * theHandler = new G4ExcitationHandler;
        theHandler->SetEvaporation(theEvaporation);
        theHandler->SetFermiModel(theFermiBreakUp);
        theHandler->SetMultiFragmentation(theMF);
        theHandler->SetMaxAandZForFermiBreakUp(12, 6);
        theHandler->SetMinEForMultiFrag(10*MeV);
	
    // Pre equilibrium stage (dummy for now)
    G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel;           
        thePreEquilib->SetExcitationHandler(theHandler);

    // make a nucleus
    G4Fancy3DNucleus * the3DNucleus = new G4Fancy3DNucleus;
    
    // make the final states
    G4ParticleScatterer * theParticleScatterer = new G4ParticleScatterer;
    // annihilation
    G4VAnnihilator * theAnihilator = new G4Annihilator;
    // elastic
    G4VElasticScatterer * theElasticScatterer = new G4ElasticScatterer;
    // inelastic
    G4VInElasticScatterer *  theInElasticScatterer = new G4InElasticScatterer;
    // cross-sections
    G4VCrossSectionBase * theCrossSectionBase = new G4CrossSectionBase;
        theParticleScatterer->SetElasticScatter(theElasticScatterer);
        theParticleScatterer->SetInElasticScatter(theInElasticScatterer);
        theParticleScatterer->SetAnnihilator(theAnihilator);
        theParticleScatterer->SetCrossSectionBase(theCrossSectionBase);
    
    // make the field propagation
    G4FieldPropagation * theFieldPropagation = new G4RKFieldIntegrator;
    
    // Here comes the kinetic model. A stripped down and highly approximate 
    // time-like cascade for the moment, to be evolved into 
    // a QMD model later, but all the bastract hooks are there
    G4HadronKineticModel * theCascade = new G4HadronKineticModel;
        theCascade->SetDeExcitation(thePreEquilib);  
	theCascade->Set3DNucleus(the3DNucleus);
	theCascade->SetFieldPropagation(theFieldPropagation);
        theCascade->SetParticleScatterer(theParticleScatterer);
	
    // here come the high energy parts
    // the string model; still not quite according to design - Explicite use of the forseen interfaces 
    // will be tested and documented in this program by beta-02 at latest.
    G4PartonStringModel * theStringModel = new G4FTFModel;
    
    theTheoModel->SetTransport(theCascade);
    theTheoModel->SetHighEnergyGenerator(theStringModel);
    theTheoModel->SetMinEnergy(20*GeV);
    theTheoModel->SetMaxEnergy(100*TeV);

//        G4DummyFragmentationModel * theFragmentation = new G4DummyFragmentationModel;
//      theStringModel->SetStringFragmentationModel(theFragmentation);
//        G4DummyGenerator * theGenerator = new G4DummyGenerator;
//      theStringModel->SetGenerator(theGenerator);
//      theStringModel->Set3DNucleus(the3DNucleus);

//    G4TheoFSGenerator * theTheoModel = new G4TheoFSGenerator;
//      G4HadronKineticModel * theCascade = new G4HadronKineticModel;
//        G4ExcitationHandler * theHandler = new G4ExcitationHandler;
//          G4StatEvaporation * theEvaporation = new G4StatEvaporation;
//            G4StatCompetitiveFission * theCompFission = new G4StatCompetitiveFission;
//            theEvaporation->SetFission(theCompFission);
//        theHandler->SetEvaporation(theEvaporation);
//          G4StatFermiBreakUp * theFermiBreakUp = new G4StatFermiBreakUp;
//        theHandler->SetFermiModel(theFermiBreakUp);
//          G4DummyMF * theMF = new G4DummyMF;
//        theHandler->SetMultiFragmentation(theMF);
//        theHandler->SetMaxAAndZForFermiBreakUp(12, 6);
//        theHandler->SetMinExEnergyForMultiFragmentation(10*MeV);
//      G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel;
//      thePreEquilib->SetExcitationHandler(theHandler);
//      theCascade->SetDeExcitation(thePreEquilib);  
//        G4Fancy3DNucleus * the3DNucleus = new G4Fancy3DNucleus;
//      theCascade->Set3DNucleus(the3DNucleus);
//        G4CrossSectionBase * theXsec = new G4CrossSectionBase;
//      theCascade->SetCrossSectionBase(theXsec);
//        G4LEProtonInelastic * theInteraction = new G4LEProtonInelastic;
//      theCascade->SetHadronicInteraction(theInteraction);
//    theTheoModel->SetTransport(theCascade);
//      G4StringModel * theStringModel = new G4StringModel;
//        G4DummyFragmentationModel * theFragmentation = new G4DummyFragmentationModel;
//      theStringModel->SetStringFragmentationModel(theFragmentation);
//        G4DummyGenerator * theGenerator = new G4DummyGenerator;
//      theStringModel->SetGenerator(theGenerator);
//      theStringModel->Set3DNucleus(the3DNucleus);
//    theTheoModel->SetHighEnergyGenerator(theStringModel);
//    theTheoModel->SetMinEnergy(20*GeV);
//    theTheoModel->SetMaxEnergy(100*TeV);
    
    G4HadronInelasticProcess *theProcesses[25];
    G4ProcessManager *theProtonProcessManager =
      theProton->GetProcessManager();
    G4ProtonInelasticProcess *theProtonInelasticProcess =
      new G4ProtonInelasticProcess();
    theProtonInelasticProcess->RegisterMe( theTheoModel );
    theProtonProcessManager->AddDiscreteProcess( theProtonInelasticProcess );
    theProcesses[0] = theProtonInelasticProcess;
   
    G4ProcessManager *theAntiProtonProcessManager =
      theAntiProton->GetProcessManager();
    G4AntiProtonInelasticProcess *theAntiProtonInelasticProcess =
      new G4AntiProtonInelasticProcess(); 
    theAntiProtonInelasticProcess->RegisterMe( theTheoModel );
    theAntiProtonProcessManager->AddDiscreteProcess( theAntiProtonInelasticProcess );
    theProcesses[1] = theAntiProtonInelasticProcess;
    
    G4ProcessManager *theNeutronProcessManager =
      theNeutron->GetProcessManager();
    G4NeutronInelasticProcess *theNeutronInelasticProcess =
      new G4NeutronInelasticProcess(); 
    theNeutronInelasticProcess->RegisterMe( theTheoModel );
    theNeutronProcessManager->AddDiscreteProcess( theNeutronInelasticProcess );
    theProcesses[2] = theNeutronInelasticProcess;
    
    G4ProcessManager *theAntiNeutronProcessManager =
      theAntiNeutron->GetProcessManager();
    G4AntiNeutronInelasticProcess *theAntiNeutronInelasticProcess =
      new G4AntiNeutronInelasticProcess();
    theAntiNeutronInelasticProcess->RegisterMe( theTheoModel );
    theAntiNeutronProcessManager->AddDiscreteProcess( theAntiNeutronInelasticProcess );
    theProcesses[3] = theAntiNeutronInelasticProcess;

    G4ProcessManager *thePionPlusProcessManager =
      thePionPlus->GetProcessManager();
    G4PionPlusInelasticProcess *thePionPlusInelasticProcess =
      new G4PionPlusInelasticProcess();
    thePionPlusInelasticProcess->RegisterMe( theTheoModel );
    thePionPlusProcessManager->AddDiscreteProcess( thePionPlusInelasticProcess );
    theProcesses[4] = thePionPlusInelasticProcess;
    
    G4ProcessManager *thePionMinusProcessManager =
      thePionMinus->GetProcessManager();
    G4PionMinusInelasticProcess *thePionMinusInelasticProcess =
      new G4PionMinusInelasticProcess();
    thePionMinusInelasticProcess->RegisterMe( theTheoModel );
    thePionMinusProcessManager->AddDiscreteProcess( thePionMinusInelasticProcess );
    theProcesses[5] = thePionMinusInelasticProcess;
    
    G4ProcessManager *theDeuteronProcessManager =
      theDeuteron->GetProcessManager();
    G4DeuteronInelasticProcess *theDeuteronInelasticProcess =
      new G4DeuteronInelasticProcess();
    theDeuteronInelasticProcess->RegisterMe( theTheoModel );
    theDeuteronProcessManager->AddDiscreteProcess( theDeuteronInelasticProcess );
    theProcesses[6] = theDeuteronInelasticProcess;
    
    G4ProcessManager *theTritonProcessManager =
      theTriton->GetProcessManager();
    G4TritonInelasticProcess *theTritonInelasticProcess =
      new G4TritonInelasticProcess();
    theTritonInelasticProcess->RegisterMe( theTheoModel );
    theTritonProcessManager->AddDiscreteProcess( theTritonInelasticProcess );
    theProcesses[7] = theTritonInelasticProcess;
    
    G4ProcessManager *theAlphaProcessManager =
      theAlpha->GetProcessManager();
    G4AlphaInelasticProcess *theAlphaInelasticProcess =
      new G4AlphaInelasticProcess();
    theAlphaInelasticProcess->RegisterMe( theTheoModel );
    theAlphaProcessManager->AddDiscreteProcess( theAlphaInelasticProcess );
    theProcesses[8] = theAlphaInelasticProcess;
    
    G4ProcessManager *theKaonPlusProcessManager =
      theKaonPlus->GetProcessManager();
    G4KaonPlusInelasticProcess *theKaonPlusInelasticProcess =
      new G4KaonPlusInelasticProcess();
    theKaonPlusInelasticProcess->RegisterMe( theTheoModel );
    theKaonPlusProcessManager->AddDiscreteProcess( theKaonPlusInelasticProcess );
    theProcesses[9] = theKaonPlusInelasticProcess;
    
    G4ProcessManager *theKaonMinusProcessManager =
      theKaonMinus->GetProcessManager();
    G4KaonMinusInelasticProcess *theKaonMinusInelasticProcess =
      new G4KaonMinusInelasticProcess();
    theKaonMinusInelasticProcess->RegisterMe( theTheoModel );
    theKaonMinusProcessManager->AddDiscreteProcess( theKaonMinusInelasticProcess );
    theProcesses[10] = theKaonMinusInelasticProcess;

    G4ProcessManager *theKaonLongProcessManager =
      theKaonLong->GetProcessManager();
    G4KaonZeroLInelasticProcess *theKaonLongInelasticProcess =
      new G4KaonZeroLInelasticProcess();
    theKaonLongInelasticProcess->RegisterMe( theTheoModel );
    theKaonLongProcessManager->AddDiscreteProcess( theKaonLongInelasticProcess );
    theProcesses[11] = theKaonLongInelasticProcess;
    
    G4ProcessManager *theKaonShortProcessManager =
      theKaonShort->GetProcessManager();
    G4KaonZeroSInelasticProcess *theKaonShortInelasticProcess =
      new G4KaonZeroSInelasticProcess();
    theKaonShortInelasticProcess->RegisterMe( theTheoModel );
    theKaonShortProcessManager->AddDiscreteProcess( theKaonShortInelasticProcess );
    theProcesses[12] = theKaonShortInelasticProcess;
    
    G4ProcessManager *theLambdaProcessManager =
      theLambda->GetProcessManager();
    G4LambdaInelasticProcess *theLambdaInelasticProcess =
      new G4LambdaInelasticProcess();
    theLambdaInelasticProcess->RegisterMe( theTheoModel );
    theLambdaProcessManager->AddDiscreteProcess( theLambdaInelasticProcess );
    theProcesses[13] = theLambdaInelasticProcess;
    
    G4ProcessManager *theAntiLambdaProcessManager =
      theAntiLambda->GetProcessManager();
    G4AntiLambdaInelasticProcess *theAntiLambdaInelasticProcess =
      new G4AntiLambdaInelasticProcess();
    theAntiLambdaInelasticProcess->RegisterMe( theTheoModel );
    theAntiLambdaProcessManager->AddDiscreteProcess( theAntiLambdaInelasticProcess );
    theProcesses[14] = theAntiLambdaInelasticProcess;
    
    G4ProcessManager *theSigmaPlusProcessManager =
      theSigmaPlus->GetProcessManager();
    G4SigmaPlusInelasticProcess *theSigmaPlusInelasticProcess =
       new G4SigmaPlusInelasticProcess();
    theSigmaPlusInelasticProcess->RegisterMe( theTheoModel );
    theSigmaPlusProcessManager->AddDiscreteProcess( theSigmaPlusInelasticProcess );
    theProcesses[15] = theSigmaPlusInelasticProcess;
    
    G4ProcessManager *theAntiSigmaPlusProcessManager =
      theAntiSigmaPlus->GetProcessManager();
    G4AntiSigmaPlusInelasticProcess *theAntiSigmaPlusInelasticProcess =
      new G4AntiSigmaPlusInelasticProcess();
    theAntiSigmaPlusInelasticProcess->RegisterMe( theTheoModel );
    theAntiSigmaPlusProcessManager->AddDiscreteProcess( theAntiSigmaPlusInelasticProcess );
    theProcesses[16] = theAntiSigmaPlusInelasticProcess;
    
    G4ProcessManager *theSigmaMinusProcessManager =
      theSigmaMinus->GetProcessManager();
    G4SigmaMinusInelasticProcess *theSigmaMinusInelasticProcess =
      new G4SigmaMinusInelasticProcess();
    theSigmaMinusInelasticProcess->RegisterMe( theTheoModel );
    theSigmaMinusProcessManager->AddDiscreteProcess( theSigmaMinusInelasticProcess );
    theProcesses[17] = theSigmaMinusInelasticProcess;
    
    G4ProcessManager *theAntiSigmaMinusProcessManager =
      theAntiSigmaMinus->GetProcessManager();
    G4AntiSigmaMinusInelasticProcess *theAntiSigmaMinusInelasticProcess =
      new G4AntiSigmaMinusInelasticProcess();
    theAntiSigmaMinusInelasticProcess->RegisterMe( theTheoModel );
    theAntiSigmaMinusProcessManager->AddDiscreteProcess( theAntiSigmaMinusInelasticProcess );
    theProcesses[18] = theAntiSigmaMinusInelasticProcess;
    
    G4ProcessManager *theXiMinusProcessManager =
      theXiMinus->GetProcessManager();
    G4XiMinusInelasticProcess *theXiMinusInelasticProcess =
      new G4XiMinusInelasticProcess();
    theXiMinusInelasticProcess->RegisterMe( theTheoModel );
    theXiMinusProcessManager->AddDiscreteProcess( theXiMinusInelasticProcess );
    theProcesses[19] = theXiMinusInelasticProcess;
    
    G4ProcessManager *theAntiXiMinusProcessManager =
      theAntiXiMinus->GetProcessManager();
    G4AntiXiMinusInelasticProcess *theAntiXiMinusInelasticProcess =
      new G4AntiXiMinusInelasticProcess();
    theAntiXiMinusInelasticProcess->RegisterMe( theTheoModel );
    theAntiXiMinusProcessManager->AddDiscreteProcess( theAntiXiMinusInelasticProcess );
    theProcesses[20] = theAntiXiMinusInelasticProcess;
    
    G4ProcessManager *theXiZeroProcessManager =
      theXiZero->GetProcessManager();
    G4XiZeroInelasticProcess *theXiZeroInelasticProcess =
      new G4XiZeroInelasticProcess();
    theXiZeroInelasticProcess->RegisterMe( theTheoModel );
    theXiZeroProcessManager->AddDiscreteProcess( theXiZeroInelasticProcess );
    theProcesses[21] = theXiZeroInelasticProcess;
    
    G4ProcessManager *theAntiXiZeroProcessManager =
      theAntiXiZero->GetProcessManager();
    G4AntiXiZeroInelasticProcess *theAntiXiZeroInelasticProcess =
      new G4AntiXiZeroInelasticProcess();
    theAntiXiZeroInelasticProcess->RegisterMe( theTheoModel );
    theAntiXiZeroProcessManager->AddDiscreteProcess( theAntiXiZeroInelasticProcess );
    theProcesses[22] = theAntiXiZeroInelasticProcess;
    
    G4ProcessManager *theOmegaMinusProcessManager =
      theOmegaMinus->GetProcessManager();
    G4OmegaMinusInelasticProcess *theOmegaMinusInelasticProcess =
      new G4OmegaMinusInelasticProcess();
    theOmegaMinusInelasticProcess->RegisterMe( theTheoModel );
    theOmegaMinusProcessManager->AddDiscreteProcess( theOmegaMinusInelasticProcess );
    theProcesses[23] = theOmegaMinusInelasticProcess;
    
    G4ProcessManager *theAntiOmegaMinusProcessManager =
      theAntiOmegaMinus->GetProcessManager();
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
    
    G4Step aStep;
    G4double meanFreePath;
    G4double incomingEnergy;
    G4int k, i, l;
    G4ParticleDefinition *currentDefinition;
    G4Track *currentTrack;
    
    G4cout << " 0) Copper      1) Lead         2) Iron      3) Tungsten" << G4endl;
    G4cout << " 4) LArgon      5) PolyStyrene  6) PbWO4     7) Oxygen" << G4endl;
    G4cout << " 8) Beryllium   9) Aluminium   10) Uranium  11) BGO" << G4endl;
    G4cout << "12) NaI        13) CsI         14) Kapton" << G4endl;
    G4cout << "Please enter the material code >> ";
    G4cin >> k;
    G4cout << G4endl;
    
    G4cout << " 0) proton        1) anti_proton    2) neutron       3) anti_neutron" << G4endl;
    G4cout << " 4) pi+           5) pi-            6) deuteron      7) triton" << G4endl;
    G4cout << " 8) alpha         9) kaon+         10) kaon-        11) kaon0L" << G4endl;
    G4cout << "12) kaon0S       13) lambda        14) anti_lambda  15) sigma+" << G4endl;
    G4cout << "16) anti_sigma+  17) sigma-        18) anti_sigma-  19) xi-" << G4endl;
    G4cout << "20) anti_xi-     21) xi0           22) anti_xi0     23) omega-" << G4endl;
    G4cout << "24) anti_omega-" << G4endl;
    G4cout << "Please enter the particle type code >> ";
    G4cin >> i;
    G4cout << G4endl;
    
    // energies of interest: 0.5, 1.0, 3.0, 12.0, 20.0 GeV

    G4cout << "Please enter the initial kinetic energy (GeV) >> ";
    G4cin >> incomingEnergy;
    G4cout << G4endl;
    incomingEnergy *= GeV/MeV;
    
    G4int nEvents;
    G4cout << "Please enter the number of events >> ";
    G4cin >> nEvents;
    G4cout << G4endl;
    
    theParticles[i]->SetCuts( 1.0*mm );
    LogicalFrame->SetMaterial( theMaterials[k] ); 
    
    // --------- Test the PostStepDoIt now  --------------
    
    G4ParticleChange *aFinalState;
    LogicalFrame->SetMaterial( theMaterials[k] ); 
    G4DynamicParticle aParticle( theParticles[i], theDirection, incomingEnergy );
    
    fprintf( g4data, "%s\n", theParticles[i]->GetParticleName() );
    fprintf( g4data, "%s\n", theMaterials[k]->GetName() );
    fprintf( g4data, "%6.1f\n", incomingEnergy/GeV );
    
    G4double currentPx, currentPy, currentPz;
    G4double kineticEnergy, totalEnergy, currentMass;
    
    G4int eventCounter = 0;
    for( l=0; l<nEvents; ++l )   // event loop
    {
      aParticle.SetDefinition( theParticles[i] );
      aParticle.SetMomentumDirection( theDirection );
      aParticle.SetKineticEnergy( incomingEnergy );
      
      G4Track *aTrack = new G4Track( &aParticle, aTime, aPosition );
      aTrack->SetTouchable( theTouchable );
      aTrack->SetStep( &aStep );
      aStep.SetTrack( aTrack );
      G4cout << "STARTING WITH EVENT NUMBER " << l+1 <<" <==================================="<<G4endl;
      aFinalState = (G4ParticleChange * ) (theProcesses[i]->PostStepDoIt( *aTrack, aStep ) );
      G4int nSec = aFinalState->GetNumberOfSecondaries();
      //delete aTrack;
      
      // prepare the analysis current quantities
      
      G4int istart, ipar;        
      if( aFinalState->GetStatusChange() == fAlive ) // current particle is still alive
      {
        istart = 0;
        G4cout <<" aParticle.GetDefinition 1"<<G4endl;
        currentDefinition = aParticle.GetDefinition();
        kineticEnergy = aFinalState->GetEnergyChange()/GeV;
        currentMass = currentDefinition->GetPDGMass()/GeV;
        totalEnergy = kineticEnergy+currentMass;
        const G4ParticleMomentum *mom = aFinalState->GetMomentumChange();
        G4double p = sqrt( kineticEnergy*kineticEnergy +
                           2.0*kineticEnergy*currentMass );
        currentPx = mom->x()*p;
        currentPy = mom->y()*p;
        currentPz = mom->z()*p;
        ipar = NameToGeant3Number( currentDefinition );
        fprintf( g4data, "%13.6E %13.6E %13.6E %13.6E %13.6E %2d %3d %6d\n",
                 kineticEnergy, totalEnergy, currentPx, currentPy, currentPz,
                 ipar, nSec+1, l );
      } else if (nSec !=0){
        istart = 1;
        currentTrack = aFinalState->GetSecondary(0);
        G4cout <<" aParticle.GetDefinition 2"<<G4endl;
        currentDefinition = currentTrack->GetDefinition();
        kineticEnergy = currentTrack->GetKineticEnergy()/GeV;
        currentMass = currentDefinition->GetPDGMass()/GeV;
        totalEnergy = kineticEnergy+currentMass;
        currentPx = currentTrack->GetMomentum().x()/GeV;
        currentPy = currentTrack->GetMomentum().y()/GeV;
        currentPz = currentTrack->GetMomentum().z()/GeV;
        ipar = NameToGeant3Number( currentDefinition );
        fprintf( g4data, "%13.6E %13.6E %13.6E %13.6E %13.6E %2d %3d %6d\n",
                 kineticEnergy, totalEnergy, currentPx, currentPy, currentPz,
                 ipar, nSec, l );
      }
      for( G4int ii=istart; ii<nSec; ++ii )
      {
        G4Track *currentTrack; 
        currentTrack = aFinalState->GetSecondary(ii);
        G4cout <<" aParticle.GetDefinition 3"<<G4endl;
        currentDefinition = currentTrack->GetDefinition();
        kineticEnergy = currentTrack->GetKineticEnergy()/GeV;
        currentMass = currentDefinition->GetPDGMass()/GeV;
        totalEnergy = kineticEnergy+currentMass;
        currentPx = currentTrack->GetMomentum().x();
        currentPy = currentTrack->GetMomentum().y();
        currentPz = currentTrack->GetMomentum().z();
        ipar = NameToGeant3Number( currentDefinition );
        fprintf( g4data, "%13.6E %13.6E %13.6E %13.6E %13.6E %2d %3d %6d\n",
                 kineticEnergy, totalEnergy, currentPx, currentPy, currentPz,
                 ipar, 0, l );
      }
      aFinalState->Clear();
    }  // end of event loop
    fclose( g4data );
    return EXIT_SUCCESS;
  }
 /* end of code */
 
