// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: CHIPS_StringModel_Analysis.cc,v 1.2 2000-09-13 07:09:05 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Johannes Peter Wellisch, 22.Apr 1997: full test-suite coded.    
// Viktor Krylov, 6.Dec 1998: some corrections.
    
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
 
#include "G4TheoFSGenerator.hh"
#include "G4StringChipsInterface.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4StringModel.hh"
#include "G4VPreCompoundModel.hh"
#include "G4FTFModel.hh"
#include "G4QGSModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

#include "G4HEProtonInelastic.hh"

#include "G4ExcitationHandler.hh"
#include "G4Evaporation.hh"
#include "G4CompetitiveFission.hh"
#include "G4FermiBreakUp.hh"
#include "G4StatMF.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4PreCompoundModel.hh"

#include "G4StringInfoDump.hh"

#include "G4DynamicParticle.hh"
#include "G4Proton.hh"
#include "G4DymmyINC.hh"

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
    G4Material *theH = new G4Material(name="Hydrogen", density=1.53*g/cm3, nEl=1); 
    theH->AddElement( elH, 1 );
    G4Material *theC = new G4Material(name="Carbon", density=1.032*g/cm3, nEl=1);
    theC->AddElement( elC, 1 );
    G4Material *theLi6 = new G4Material(name="Li6", density=2.70*g/cm3, nEl=1);
    G4int nIso;
    G4Element *elLi6  = new G4Element(name="Li6", symbol="Li", nIso = 1);
    G4Isotope * isoLi6 = new G4Isotope(name="Li6", 3, 6, a=6.*g/mole);
    elLi6->AddIsotope(isoLi6, 1);
    theLi6->AddElement(elLi6 , 1 );
    G4Material *theTa = new G4Material(name="Tantalum", density=1.53*g/cm3, nEl=1); 
    G4Element *elTa = new G4Element(name="Tantalum", symbol="Ta", iz=73., a=181*g/mole);
    theTa->AddElement( elTa, 1 );
    
    G4int numberOfMaterials = 19;
    G4Material *theMaterials[19];
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
    theMaterials[15] = theH;
    theMaterials[16] = theC;
    theMaterials[17] = theLi6;
    theMaterials[18] = theTa;
   
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

    // This is a parameterized model
    G4HEProtonInelastic * theParaModel = new G4HEProtonInelastic;

    // This is the preparation for pre-compound
    G4Evaporation * theEvaporation = new G4Evaporation;
    G4FermiBreakUp * theFermiBreakUp = new G4FermiBreakUp;
    G4StatMF * theMF = new G4StatMF;

    G4ExcitationHandler * theHandler = new G4ExcitationHandler;
        theHandler->SetEvaporation(theEvaporation);
        theHandler->SetFermiModel(theFermiBreakUp);
        theHandler->SetMultiFragmentation(theMF);
        theHandler->SetMaxAandZForFermiBreakUp(12, 6);
        theHandler->SetMinEForMultiFrag(3*MeV);
	
    G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel(theHandler);

    G4GeneratorPrecompoundInterface * theCascade = new G4GeneratorPrecompoundInterface;
            theCascade->SetDeExcitation(thePreEquilib);  
   
    // this will be the model class for high energies
    G4TheoFSGenerator * theTheoModel = new G4TheoFSGenerator;
       
    // all models for treatment of thermal nucleus 

    // Evaporation logic
	
    // Pre equilibrium stage (dummy for now)

    // make a nucleus
    G4Fancy3DNucleus * the3DNucleus = new G4Fancy3DNucleus;
    
    // here come the high energy parts
    // the string model; still not quite according to design - Explicite use of the forseen interfaces 
    // will be tested and documented in this program by beta-02 at latest.
    G4VPartonStringModel * theStringModel;
    G4int choice;
    do {
      G4cout << "Please select String Model to test(1-FTF, 2-QGS): " << G4std::flush;
      G4cin >> choice;
      if(choice==1) theStringModel = new G4FTFModel;
      if(choice==2) theStringModel = new G4QGSModel;
    }while(choice!=1&&choice!=2);

    theTheoModel->SetHighEnergyGenerator(theStringModel);
    theTheoModel->SetMinEnergy(20*GeV);
    theTheoModel->SetMaxEnergy(100*TeV);

      G4cout << "Please select fragmentation model [1=lund, 2=qgsm] "<<G4endl;
      G4int ifragmodel;
      G4cin >> ifragmodel;
      G4VLongitudinalStringDecay * theFragmentation;
      if(ifragmodel == 1)
        theFragmentation = new G4LundStringFragmentation;
      else if(ifragmodel == 2)
        theFragmentation = new G4QGSMFragmentation;
      else G4Exception("G4StringModelTest: unknown fragmentation model");
      G4ExcitedStringDecay theStringDecay(theFragmentation);
      theStringModel->SetFragmentationModel(&theStringDecay);
    
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
      G4cout << " 0) Copper      1) Lead         2) Iron      3) Tungsten" << G4endl;
      G4cout << " 4) LArgon      5) PolyStyrene  6) PbWO4     7) Oxygen" << G4endl;
      G4cout << " 8) Beryllium   9) Aluminium   10) Uranium  11) BGO" << G4endl;
      G4cout << "12) NaI        13) CsI         14) Kapton   15) Hydrogen" << G4endl;
      G4cout << "16) Carbon     17) 6_3_Lithium 18) Tantalum "<< G4endl;
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

    // energies of interest: 0.5, 1.0, 3.0, 12.0, 20.0 GeV

    G4cout << "Please enter the initial kinetic energy (GeV): "<< G4std::flush;
    G4cin >> incomingEnergy;
    incomingEnergy *= GeV/MeV;

    theTheoModel->SetTransport(theCascade);
    G4StringChipsInterface * theChipsCascade = new G4StringChipsInterface();
    theTheoModel->SetTransport(theChipsCascade);
//    G4StringInfoDump * theDummyCascade = new G4StringInfoDump;
//    G4DymmyINC * theDummyCascade = new G4DymmyINC;
//    theTheoModel->SetTransport(theDummyCascade);
    
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
            delete aSec;
	    delete second;
          }
          delete aTrack;
	  aFinalState->Clear();
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
 
