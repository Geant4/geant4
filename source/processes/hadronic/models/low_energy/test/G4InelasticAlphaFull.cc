// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4InelasticAlphaFull.cc,v 1.2 1999-12-15 14:53:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Johannes Peter Wellisch, 22.Apr 1997: full test-suite coded.    
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"
 
#include "G4Material.hh"
 
#include "G4GRSVolume.hh"
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
 
#include "G4LEAlphaInelastic.hh"
#include "G4LEAntiLambdaInelastic.hh"
#include "G4LEAntiNeutronInelastic.hh"
#include "G4LEAntiOmegaMinusInelastic.hh"
#include "G4LEAntiProtonInelastic.hh"
#include "G4LEAntiSigmaMinusInelastic.hh"
#include "G4LEAntiSigmaPlusInelastic.hh"
#include "G4LEAntiXiMinusInelastic.hh"
#include "G4LEAntiXiZeroInelastic.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4LEKaonMinusInelastic.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4LEKaonZeroLInelastic.hh"
#include "G4LEKaonZeroSInelastic.hh"
#include "G4LELambdaInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4LEOmegaMinusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LESigmaMinusInelastic.hh"
#include "G4LESigmaPlusInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEXiMinusInelastic.hh"
#include "G4LEXiZeroInelastic.hh"

#include "G4DynamicParticle.hh"
#include "G4Proton.hh"
#include "G4ProcessManager.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"

#include "G4Step.hh"

#include "g4templates.hh"
 
 // forward declarations
 
 G4int sortEnergies( const double Px, const double Py, const double Pz,
                     const double Ekin, double* sortedPx, double* sortedPy,
                     double* sortedPz, double* sortedE );
 
 // here comes the code
 
 G4int sortEnergies( const double Px, const double Py, const double Pz,
                     const double Ekin, double* sortedPx, double* sortedPy,
                     double* sortedPz, double* sortedE)
  {
    for( int i=0; i<10; ++i ) {
      if( abs(Ekin) > sortedE[i] ) {
        sortedE[i]  = Ekin;
        sortedPx[i] = Px;
        sortedPy[i] = Py;
        sortedPz[i] = Pz;
        return 1;
      }
    }
    return 0;
  }
 
 int main()
  {
    G4cout.setf( G4std::ios::scientific, G4std::ios::floatfield );
    G4std::ofstream outFile( "InelasticAlpha.listing.GetMeanFreePath", G4std::ios::out);
    outFile.setf( G4std::ios::scientific, G4std::ios::floatfield );
    G4std::ofstream outFile1( "InelasticAlpha.listing.DoIt", G4std::ios::out);
    outFile1.setf( G4std::ios::scientific, G4std::ios::floatfield );

    G4String name, symbol;
    G4double a, iz, z, density;
    G4int nEl;
    
    G4Material* theCu = new G4Material(name="Copper", density=8.96*g/cm3, nEl=1);
    G4Element* elCu = new G4Element(name="Copper", symbol="Cu", iz=29., a=63.55*g/mole);
    theCu->AddElement(elCu, 1);
    G4Material* thePb = new G4Material(name="Lead", density=11.35*g/cm3, nEl=1);
    G4Element* elPb = new G4Element(name="Lead", symbol="Pb", iz=82., a=207.19*g/mole);
    thePb->AddElement(elPb, 1);
    G4Material* theFe = new G4Material(name="Iron", density=7.87*g/cm3, nEl=1);
    G4Element* elFe = new G4Element(name="Iron", symbol="Fe", iz=26., a=55.85*g/mole);
    theFe->AddElement(elFe, 1);
    G4Material* theW  = new G4Material(name="Tungsten", density=19.30*g/cm3, nEl=1);
    G4Element* elW  = new G4Element(name="Tungston", symbol="W", iz=74., a=183.85*g/mole);
    theW->AddElement(elW, 1);
    G4Material* theLAr= new G4Material(name="LArgon", density=1.393*g/cm3, nEl=1);
    G4Element* elAr  = new G4Element(name="Argon", symbol="Ar", iz=18., a=39.95*g/mole);
    theLAr->AddElement(elAr, 1);
    G4Material* thePS = new G4Material(name="PolyStyrene", density=1.032*g/cm3, nEl=2);
    G4Element* elC = new G4Element(name="Carbon", symbol="C", iz=6., a=12.01*g/mole);
    G4Element* elH = new G4Element(name="Hydrogen", symbol="H", iz=1., a=1.01*g/mole);
    thePS->AddElement(elC, 8);
    thePS->AddElement(elH, 8);
    G4Material* thePbWO4 = new G4Material(name="PbWO4", density=12.0*g/cm3, nEl=3);
    // approximative number
    G4Element* elO = new G4Element(name="Oxygen", symbol="O", iz=8., a=15.9994*g/mole);
    thePbWO4->AddElement(elPb, 1);
    thePbWO4->AddElement(elW,  1);
    thePbWO4->AddElement(elO,  4);
    // approximate numbers for O
    G4Material* theO = new G4Material(name="Oxygen", density=1.1*g/cm3, nEl=1);
    theO->AddElement(elO,  1);
    G4Material* theBe = new G4Material(name="Beryllium", density=1.848*g/cm3, nEl=1);
    G4Element* elBe  = new G4Element(name="Beryllium", symbol="Be", iz=4., a=9.01*g/mole);
    theBe->AddElement(elBe, 1);
    G4Material* theAl = new G4Material(name="Aluminium", density=2.70*g/cm3, nEl=1);
    G4Element* elAl  = new G4Element(name="Aluminium", symbol="Al", iz=13., a=26.98*g/mole);
    theAl->AddElement(elAl, 1);
    G4Material* theU = new G4Material(name="Uranium", density=18.95*g/cm3, nEl=1);
    G4Element* elU  = new G4Element(name="Uranium", symbol="U", iz=92., a=238.03*g/mole);
    theU->AddElement(elU, 1);
    G4Material* theBGO = new G4Material(name="BGO", density=2.15*g/cm3, nEl=3);
    G4Element* elBi = new G4Element(name="Bismuth", symbol="Bi", iz=83., a=208.98*g/mole);
    G4Element* elGe = new G4Element(name="Germanium", symbol="Ge", iz=32., a=72.59*g/mole);
    theBGO->AddElement(elBi, 4);
    theBGO->AddElement(elGe, 3);
    theBGO->AddElement(elO, 12);
    G4Material* theNaI = new G4Material(name="NaI", density=3.67*g/cm3, nEl=2);
    G4Element* elNa = new G4Element(name="Sodium", symbol="Na", iz=11., a=22.990*g/mole);
    G4Element* elI = new G4Element(name="Iodine", symbol="I", iz=53., a=126.904*g/mole);
    theNaI->AddElement(elNa, 1);
    theNaI->AddElement(elI, 1);
    G4Material* theCsI = new G4Material(name="CsI", density=4.53*g/cm3, nEl=2);
    G4Element* elCs = new G4Element(name="Cesium", symbol="Cs", iz=55., a=132.905*g/mole);
    theCsI->AddElement(elCs, 1);
    theCsI->AddElement(elI, 1);
    G4Material* theKapton = new G4Material(name="Kapton", density=1.53*g/cm3, nEl=4); 
    // formula: private communications, see mail.
    theKapton->AddElement(elC, 22);
    theKapton->AddElement(elH, 10);
    theKapton->AddElement(elO, 5);
    G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a=14.007*g/mole);
    theKapton->AddElement(elN, 2);
    
    G4int numberOfMaterials=15;
    G4Material* theMaterials[15];
    theMaterials[0]=theCu;
    theMaterials[1]=thePb;
    theMaterials[2]=theFe;
    theMaterials[3]=theW;
    theMaterials[4]=theLAr;
    theMaterials[5]=thePS;
    theMaterials[6]=thePbWO4;
    theMaterials[7]=theO;
    theMaterials[8]=theBe;
    theMaterials[9]=theAl;
    theMaterials[10]=theU;
    theMaterials[11]=theBGO;
    theMaterials[12]=theNaI;
    theMaterials[13]=theCsI;
    theMaterials[14]=theKapton;
    
    // ----------- here all material have been defined -----------
    
//    G4Element::DumpInfo(); 
//    G4Material::DumpInfo();
    
    // ----------- the following is needed for building a track...... ------------
    
    static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
    G4int imat = 0;   
    G4Box* theFrame = new G4Box ("Frame",10*m, 10*m, 10*m);
    
    G4LogicalVolume* LogicalFrame = new G4LogicalVolume(theFrame,
                                                        (*theMaterialTable)(imat),
                                                        "LFrame", 0, 0, 0);
    
    G4PVPlacement* PhysicalFrame = new G4PVPlacement(0,G4ThreeVector(),
                                                     "PFrame",LogicalFrame,0,false,0);
    G4RotationMatrix theNull;
    G4ThreeVector theCenter(0,0,0);
    G4GRSVolume * theTouchable = new G4GRSVolume(PhysicalFrame, &theNull, theCenter);
    // ----------- now get all particles of interest ---------
    G4int numberOfParticles = 25;
   G4ParticleDefinition* theParticles[25];
   G4ParticleDefinition* theProton = G4Proton::ProtonDefinition();
   theParticles[0]=theProton;
   G4ParticleDefinition* theAntiProton = G4AntiProton::AntiProtonDefinition();
   theParticles[1]=theAntiProton;
   G4ParticleDefinition* theNeutron = G4Neutron::NeutronDefinition();
   theParticles[2]=theNeutron;
   G4ParticleDefinition* theAntiNeutron = G4AntiNeutron::AntiNeutronDefinition();
   theParticles[3]=theAntiNeutron;
   G4ParticleDefinition* thePionPlus = G4PionPlus::PionPlusDefinition();
   theParticles[4]=thePionPlus;
   G4ParticleDefinition* thePionMinus = G4PionMinus::PionMinusDefinition();
   theParticles[5]=thePionMinus;
   G4ParticleDefinition* theDeuteron = G4Deuteron::DeuteronDefinition();
   theParticles[6]=theDeuteron;
   G4ParticleDefinition* theTriton = G4Triton::TritonDefinition();
   theParticles[7]=theTriton;
   G4ParticleDefinition* theAlpha = G4Alpha::AlphaDefinition();
   theParticles[8]=theAlpha;
   G4ParticleDefinition* theKaonPlus = G4KaonPlus::KaonPlusDefinition();
   theParticles[9]=theKaonPlus;
   G4ParticleDefinition* theKaonMinus = G4KaonMinus::KaonMinusDefinition();
   theParticles[10]=theKaonMinus;
   G4ParticleDefinition* theKaonLong = G4KaonZeroLong::KaonZeroLongDefinition();
   theParticles[11]=theKaonLong;
   G4ParticleDefinition* theKaonShort = G4KaonZeroShort::KaonZeroShortDefinition();
   theParticles[12]=theKaonShort;
   G4ParticleDefinition* theLambda = G4Lambda::LambdaDefinition();
   theParticles[13]=theLambda;
   G4ParticleDefinition* theAntiLambda = G4AntiLambda::AntiLambdaDefinition();
   theParticles[14]=theAntiLambda;
   G4ParticleDefinition* theSigmaPlus = G4SigmaPlus::SigmaPlusDefinition();
   theParticles[15]=theSigmaPlus;
   G4ParticleDefinition* theAntiSigmaPlus = G4AntiSigmaPlus::AntiSigmaPlusDefinition();
   theParticles[16]=theAntiSigmaPlus;
   G4ParticleDefinition* theSigmaMinus = G4SigmaMinus::SigmaMinusDefinition();
   theParticles[17]=theSigmaMinus;
   G4ParticleDefinition* theAntiSigmaMinus = G4AntiSigmaMinus::AntiSigmaMinusDefinition();
   theParticles[18]=theAntiSigmaMinus;
   G4ParticleDefinition* theXiMinus = G4XiMinus::XiMinusDefinition();
   theParticles[19]=theXiMinus;
   G4ParticleDefinition* theAntiXiMinus = G4AntiXiMinus::AntiXiMinusDefinition();
   theParticles[20]=theAntiXiMinus;
   G4ParticleDefinition* theXiZero = G4XiZero::XiZeroDefinition();
   theParticles[21]=theXiZero;
   G4ParticleDefinition* theAntiXiZero = G4AntiXiZero::AntiXiZeroDefinition();
   theParticles[22]=theAntiXiZero;
   G4ParticleDefinition* theOmegaMinus = G4OmegaMinus::OmegaMinusDefinition();
   theParticles[23]=theOmegaMinus;
   G4ParticleDefinition* theAntiOmegaMinus = G4AntiOmegaMinus::AntiOmegaMinusDefinition();
   theParticles[24]=theAntiOmegaMinus;
   
   //------ here all the particles are Done ----------
   G4cout << "Done with all the particles" << G4endl;
   G4cout << "Starting process definitions" << G4endl;
   //--------- Processes definitions ---------
   G4HadronInelasticProcess* theProcesses[25];
   
   G4ProcessManager* theProtonProcessManager = theProton->GetProcessManager();
   G4ProtonInelasticProcess theProtonInelasticProcess; 
   G4LEProtonInelastic theProtonInelastic;
   theProtonInelasticProcess.RegisterMe(&theProtonInelastic);
   theProtonProcessManager->AddDiscreteProcess(&theProtonInelasticProcess);
   theProcesses[0] = &theProtonInelasticProcess;
   
   G4ProcessManager* theAntiProtonProcessManager = theAntiProton->GetProcessManager();
   G4AntiProtonInelasticProcess theAntiProtonInelasticProcess; 
   G4LEAntiProtonInelastic theAntiProtonInelastic;
   theAntiProtonInelasticProcess.RegisterMe(&theAntiProtonInelastic);
   theAntiProtonProcessManager->AddDiscreteProcess(&theAntiProtonInelasticProcess);
   theProcesses[1] = &theAntiProtonInelasticProcess;
   
   G4ProcessManager* theNeutronProcessManager = theNeutron->GetProcessManager();
   G4NeutronInelasticProcess theNeutronInelasticProcess; 
   G4LENeutronInelastic theNeutronInelastic;
   theNeutronInelasticProcess.RegisterMe(&theNeutronInelastic);
   theNeutronProcessManager->AddDiscreteProcess(&theNeutronInelasticProcess);
   theProcesses[2] = &theNeutronInelasticProcess;
   
   G4ProcessManager* theAntiNeutronProcessManager = theAntiNeutron->GetProcessManager();
   G4AntiNeutronInelasticProcess theAntiNeutronInelasticProcess;
   G4LEAntiNeutronInelastic theAntiNeutronInelastic;
   theAntiNeutronInelasticProcess.RegisterMe(&theAntiNeutronInelastic);
   theAntiNeutronProcessManager->AddDiscreteProcess(&theAntiNeutronInelasticProcess);
   theProcesses[3] = &theAntiNeutronInelasticProcess;

   G4ProcessManager* thePionPlusProcessManager = thePionPlus->GetProcessManager();
   G4PionPlusInelasticProcess thePionPlusInelasticProcess;
   G4LEPionPlusInelastic thePionPlusInelastic;
   thePionPlusInelasticProcess.RegisterMe(&thePionPlusInelastic);
   thePionPlusProcessManager->AddProcess(&thePionPlusInelasticProcess,-1,-1,0);
   theProcesses[4] = &thePionPlusInelasticProcess;

   G4ProcessManager* thePionMinusProcessManager = thePionMinus->GetProcessManager();
   G4PionMinusInelasticProcess thePionMinusInelasticProcess;
   G4LEPionMinusInelastic thePionMinusInelastic;
   thePionMinusInelasticProcess.RegisterMe(&thePionMinusInelastic);
   thePionMinusProcessManager->AddProcess(&thePionMinusInelasticProcess,-1,-1,0);
   theProcesses[5] = &thePionMinusInelasticProcess;
   
   G4ProcessManager* theDeuteronProcessManager = theDeuteron->GetProcessManager();
   G4DeuteronInelasticProcess theDeuteronInelasticProcess;
   G4LEDeuteronInelastic theDeuteronInelastic;
   theDeuteronInelasticProcess.RegisterMe(&theDeuteronInelastic);
   theDeuteronProcessManager->AddProcess(&theDeuteronInelasticProcess,-1,-1,0);
   theProcesses[6] = &theDeuteronInelasticProcess;

   G4ProcessManager* theTritonProcessManager = theTriton->GetProcessManager();
   G4TritonInelasticProcess theTritonInelasticProcess;
   G4LETritonInelastic theTritonInelastic;
   theTritonInelasticProcess.RegisterMe(&theTritonInelastic);
   theTritonProcessManager->AddProcess(&theTritonInelasticProcess,-1,-1,0);
   theProcesses[7] = &theTritonInelasticProcess;

   G4ProcessManager* theAlphaProcessManager = theAlpha->GetProcessManager();
   G4AlphaInelasticProcess theAlphaInelasticProcess;
   G4LEAlphaInelastic theAlphaInelastic;
   theAlphaInelasticProcess.RegisterMe(&theAlphaInelastic);
   theAlphaProcessManager->AddProcess(&theAlphaInelasticProcess,-1,-1,0);
   theProcesses[8] = &theAlphaInelasticProcess;

   G4ProcessManager* theKaonPlusProcessManager = theKaonPlus->GetProcessManager();
   G4KaonPlusInelasticProcess theKaonPlusInelasticProcess;
   G4LEKaonPlusInelastic theKaonPlusInelastic;
   theKaonPlusInelasticProcess.RegisterMe(&theKaonPlusInelastic);
   theKaonPlusProcessManager->AddProcess(&theKaonPlusInelasticProcess,-1,-1,0);
   theProcesses[9] = &theKaonPlusInelasticProcess;

   G4ProcessManager* theKaonMinusProcessManager = theKaonMinus->GetProcessManager();
   G4KaonMinusInelasticProcess theKaonMinusInelasticProcess;
   G4LEKaonMinusInelastic theKaonMinusInelastic;
   theKaonMinusInelasticProcess.RegisterMe(&theKaonMinusInelastic);
   theKaonMinusProcessManager->AddProcess(&theKaonMinusInelasticProcess,-1,-1,0);
   theProcesses[10] = &theKaonMinusInelasticProcess;

   G4ProcessManager* theKaonLongProcessManager = theKaonLong->GetProcessManager();
   G4KaonZeroLInelasticProcess theKaonLongInelasticProcess; 
   G4LEKaonZeroLInelastic theKaonLongInelastic;
   theKaonLongInelasticProcess.RegisterMe(&theKaonLongInelastic);
   theKaonLongProcessManager->AddProcess(&theKaonLongInelasticProcess,-1,-1,0);
   theProcesses[11] = &theKaonLongInelasticProcess;
   
   G4ProcessManager* theKaonShortProcessManager = theKaonShort->GetProcessManager();
   G4KaonZeroSInelasticProcess theKaonShortInelasticProcess;
   G4LEKaonZeroSInelastic theKaonShortInelastic;
   theKaonShortInelasticProcess.RegisterMe(&theKaonShortInelastic);
   theKaonShortProcessManager->AddProcess(&theKaonShortInelasticProcess,-1,-1,0);
   theProcesses[12] = &theKaonShortInelasticProcess;
   
   G4ProcessManager* theLambdaProcessManager = theLambda->GetProcessManager();
   G4LambdaInelasticProcess theLambdaInelasticProcess;
   G4LELambdaInelastic theLambdaInelastic;
   theLambdaInelasticProcess.RegisterMe(&theLambdaInelastic);
   theLambdaProcessManager->AddProcess(&theLambdaInelasticProcess,-1,-1,0);
   theProcesses[13] = &theLambdaInelasticProcess;

   G4ProcessManager* theAntiLambdaProcessManager = theAntiLambda->GetProcessManager();
   G4AntiLambdaInelasticProcess theAntiLambdaInelasticProcess;
   G4LEAntiLambdaInelastic theAntiLambdaInelastic;
   theAntiLambdaInelasticProcess.RegisterMe(&theAntiLambdaInelastic);
   theAntiLambdaProcessManager->AddProcess(&theAntiLambdaInelasticProcess,-1,-1,0);
   theProcesses[14] = &theAntiLambdaInelasticProcess;

   G4ProcessManager* theSigmaPlusProcessManager = theSigmaPlus->GetProcessManager();
   G4SigmaPlusInelasticProcess theSigmaPlusInelasticProcess;
   G4LESigmaPlusInelastic theSigmaPlusInelastic;
   theSigmaPlusInelasticProcess.RegisterMe(&theSigmaPlusInelastic);
   theSigmaPlusProcessManager->AddProcess(&theSigmaPlusInelasticProcess,-1,-1,0);
   theProcesses[15] = &theSigmaPlusInelasticProcess;

   G4ProcessManager* theAntiSigmaPlusProcessManager = theAntiSigmaPlus->GetProcessManager();
   G4AntiSigmaPlusInelasticProcess theAntiSigmaPlusInelasticProcess;
   G4LEAntiSigmaPlusInelastic theAntiSigmaPlusInelastic;
   theAntiSigmaPlusInelasticProcess.RegisterMe(&theAntiSigmaPlusInelastic);
   theAntiSigmaPlusProcessManager->AddProcess(&theAntiSigmaPlusInelasticProcess,-1,-1,0);
   theProcesses[16] = &theAntiSigmaPlusInelasticProcess;

   G4ProcessManager* theSigmaMinusProcessManager = theSigmaMinus->GetProcessManager();
   G4SigmaMinusInelasticProcess theSigmaMinusInelasticProcess;
   G4LESigmaMinusInelastic theSigmaMinusInelastic;
   theSigmaMinusInelasticProcess.RegisterMe(&theSigmaMinusInelastic);
   theSigmaMinusProcessManager->AddProcess(&theSigmaMinusInelasticProcess,-1,-1,0);
   theProcesses[17] = &theSigmaMinusInelasticProcess;

   G4ProcessManager* theAntiSigmaMinusProcessManager = theAntiSigmaMinus->GetProcessManager();
   G4AntiSigmaMinusInelasticProcess theAntiSigmaMinusInelasticProcess;
   G4LEAntiSigmaMinusInelastic theAntiSigmaMinusInelastic;
   theAntiSigmaMinusInelasticProcess.RegisterMe(&theAntiSigmaMinusInelastic);
   theAntiSigmaMinusProcessManager->AddProcess(&theAntiSigmaMinusInelasticProcess,-1,-1,0);
   theProcesses[18] = &theAntiSigmaMinusInelasticProcess;

   G4ProcessManager* theXiMinusProcessManager = theXiMinus->GetProcessManager();
   G4XiMinusInelasticProcess theXiMinusInelasticProcess;
   G4LEXiMinusInelastic theXiMinusInelastic;
   theXiMinusInelasticProcess.RegisterMe(&theXiMinusInelastic);
   theXiMinusProcessManager->AddProcess(&theXiMinusInelasticProcess,-1,-1,0);
   theProcesses[19] = &theXiMinusInelasticProcess;

   G4ProcessManager* theAntiXiMinusProcessManager = theAntiXiMinus->GetProcessManager();
   G4AntiXiMinusInelasticProcess theAntiXiMinusInelasticProcess;
   G4LEAntiXiMinusInelastic theAntiXiMinusInelastic;
   theAntiXiMinusInelasticProcess.RegisterMe(&theAntiXiMinusInelastic);
   theAntiXiMinusProcessManager->AddProcess(&theAntiXiMinusInelasticProcess,-1,-1,0);
   theProcesses[20] = &theAntiXiMinusInelasticProcess;

   G4ProcessManager* theXiZeroProcessManager = theXiZero->GetProcessManager();
   G4XiZeroInelasticProcess theXiZeroInelasticProcess;
   G4LEXiZeroInelastic theXiZeroInelastic;
   theXiZeroInelasticProcess.RegisterMe(&theXiZeroInelastic);
   theXiZeroProcessManager->AddProcess(&theXiZeroInelasticProcess,-1,-1,0);
   theProcesses[21] = &theXiZeroInelasticProcess;

   G4ProcessManager* theAntiXiZeroProcessManager = theAntiXiZero->GetProcessManager();
   G4AntiXiZeroInelasticProcess theAntiXiZeroInelasticProcess;
   G4LEAntiXiZeroInelastic theAntiXiZeroInelastic;
   theAntiXiZeroInelasticProcess.RegisterMe(&theAntiXiZeroInelastic);
   theAntiXiZeroProcessManager->AddProcess(&theAntiXiZeroInelasticProcess,-1,-1,0);
   theProcesses[22] = &theAntiXiZeroInelasticProcess;

   G4ProcessManager* theOmegaMinusProcessManager = theOmegaMinus->GetProcessManager();
   G4OmegaMinusInelasticProcess theOmegaMinusInelasticProcess;
   G4LEOmegaMinusInelastic theOmegaMinusInelastic;
   theOmegaMinusInelasticProcess.RegisterMe(&theOmegaMinusInelastic);
   theOmegaMinusProcessManager->AddProcess(&theOmegaMinusInelasticProcess,-1,-1,0);
   theProcesses[23] = &theOmegaMinusInelasticProcess;

   G4ProcessManager* theAntiOmegaMinusProcessManager
     = theAntiOmegaMinus->GetProcessManager();
   G4AntiOmegaMinusInelasticProcess theAntiOmegaMinusInelasticProcess;
   G4LEAntiOmegaMinusInelastic theAntiOmegaMinusInelastic;
   theAntiOmegaMinusInelasticProcess.RegisterMe(&theAntiOmegaMinusInelastic);
   theAntiOmegaMinusProcessManager->AddDiscreteProcess(&theAntiOmegaMinusInelasticProcess);
   theProcesses[24] = &theAntiOmegaMinusInelasticProcess;

   G4ForceCondition* condition = new G4ForceCondition;
   *condition = NotForced;

   G4cout << "Done with all the process definitions"<<G4endl;
   //   G4cout << "Building the CrossSection Tables. This will take a while"<<G4endl;
   
   // ------- Build the CrossSection Tables ------- this design can be impoved ------
   
   //   theProton->SetCuts(1.*mm);
   //   theAntiProton->SetCuts(1.*mm);
   //   theNeutron->SetCuts(1.*mm);
   //   theAntiNeutron->SetCuts(1.*mm);
   //   thePionPlus->SetCuts(1.*mm);
   //   thePionMinus->SetCuts(1.*mm);
   //   theDeuteron->SetCuts(1.*mm);
   //   theTriton->SetCuts(1.*mm);
   //   theAlpha->SetCuts(1.*mm);
   //   theKaonPlus->SetCuts(1.*mm);
   //   theKaonMinus->SetCuts(1.*mm);
   //   theKaonLong->SetCuts(1.*mm);
   //   theKaonShort->SetCuts(1.*mm);
   //   theLambda->SetCuts(1.*mm);
   //   theAntiLambda->SetCuts(1.*mm);
   //   theSigmaPlus->SetCuts(1.*mm);
   //   theAntiSigmaPlus->SetCuts(1.*mm);
   //   theSigmaMinus->SetCuts(1.*mm);
   //   theAntiSigmaMinus->SetCuts(1.*mm);
   //   theXiMinus->SetCuts(1.*mm);
   //   theAntiXiMinus->SetCuts(1.*mm);
   //   theXiZero->SetCuts(1.*mm);
   //   theAntiXiZero->SetCuts(1.*mm);
   //   theOmegaMinus->SetCuts(1.*mm);
   //   theAntiOmegaMinus->SetCuts(1.*mm);
   
   //   G4cout << "Done with the CrossSectionTables" << G4endl;
   // ----------- define energies of interest ----------------
   
   int numberOfEnergies = 5;
   G4double theEnergies[5] = { 0.5*GeV, 1.0*GeV, 3.0*GeV, 12.0*GeV, 20.0*GeV };
   
   // ----------- needed to Build a Step, and a track -------------------
   
   G4ParticleMomentum theDirection(0.,0.,1.);
   G4ThreeVector aPosition(0.,0.,0.);
   G4double aTime = 0.0;
   
   // --------- Test the GetMeanFreePath
   
   G4Step aStep;
   G4double meanFreePath;
   G4double incomingEnergy;
   G4int k, i, l, hpw = 0;
   
   G4cout << "Test the GetMeanFreePath, and DoIt: please enter the material" << G4endl;
   G4cout << "Candidates are:" << G4endl;
   for( G4int k1=0; k1<numberOfMaterials; ++k1 )
     G4cout << k1 << ") " << theMaterials[k1]->GetName() << G4endl;
   G4cout << "Test the GetMeanFreePath, and DoIt: please enter the particle type"<<G4endl;
   G4cout << "Candidates are:"<<G4endl;
   for( G4int i1=0; i1<numberOfParticles; ++i1 )
     G4cout << i1 << ") " << theParticles[i1]->GetParticleName() << G4endl;
   G4cout << "running" << G4endl;
   theParticles[i]->SetCuts( 1.0*mm );
   for ( i=0; i<numberOfParticles; i++)
   {
     outFile << G4endl
             << "New particle type: " << theParticles[i]->GetParticleName() << G4endl;
     for ( k=0; k<numberOfMaterials; k++)
     {
       outFile << G4endl
               << "Entering Material " << theMaterials[k]->GetName() << " for particle "
               << theParticles[i]->GetParticleName() << G4endl;
       LogicalFrame->SetMaterial( theMaterials[k] ); 
       for( G4int j=0; j<numberOfEnergies; ++j ) {
         incomingEnergy = theEnergies[j];
         if( i>=6 && i<=8 )incomingEnergy /= 201.;
         G4DynamicParticle* aParticle =
           new G4DynamicParticle( theParticles[i], theDirection, incomingEnergy );
         G4Track* aTrack = new G4Track( aParticle, aTime, aPosition );
         aTrack->SetTouchable(theTouchable);
         aStep.SetTrack( aTrack );
         outFile << "  " << incomingEnergy/GeV << " GeV";
         meanFreePath = theProcesses[i]->GetMeanFreePath( *aTrack, -1., condition );
         outFile<< "            " << meanFreePath*cm << G4endl ;
         delete aParticle;
         delete aTrack;
       }  // energy loop
     }  // material loop
   }  //particle loop
   outFile << " " << G4endl;
   
   // --------- Test the PostStepDoIt now, 10 events each --------------
   G4cout << "Entering the DoIt loops!!!!!"<< G4endl;
   G4VParticleChange* aFinalState;
   G4cout << "Now test the DoIt: please enter the number of events"<<G4endl;
   G4int ll0=0;
   G4cin >> ll0;
   G4cout <<"Now debug the DoIt: enter the problem event number"<< G4endl;
   G4int debugThisOne;
   G4cin >> debugThisOne;
   for (i=0; i<numberOfParticles; i++)
   {
     outFile << G4endl
             << "New particle type: " << theParticles[i]->GetParticleName()
             << " " << i << G4endl;
     for ( G4int k=0; k<numberOfMaterials; k++)
     {
       outFile << G4endl << "Entering Material " << theMaterials[k]->GetName()
               << " for particle " << theParticles[i]->GetParticleName() << G4endl;
       LogicalFrame->SetMaterial(theMaterials[k]); 
       for( G4int j=0; j<numberOfEnergies; ++j )
       {
         for( G4int l=0; l<ll0; ++l )
         {
           incomingEnergy=theEnergies[j];
           if( i>=6 && i<=8 )incomingEnergy /= 201.;
           G4DynamicParticle* aParticle =
             new G4DynamicParticle( theParticles[i], theDirection, incomingEnergy );
           G4Track* aTrack = new G4Track( aParticle, aTime, aPosition );
           aTrack->SetTouchable(theTouchable);
           aStep.SetTrack( aTrack );
           
           G4cout << "EVENTCOUNTER=" << ++hpw
                << " current energy: " << incomingEnergy
                << " of particle " << aParticle->GetDefinition()->GetParticleName() 
                << " in material " << theMaterials[k]->GetName() << G4endl;
           if (hpw==debugThisOne)
           {
            debugThisOne+=0;
           }
           aFinalState = theProcesses[i]->PostStepDoIt( *aTrack, aStep );
           delete aParticle;
           delete aTrack;
           // -------------- now to the analysis ----------------
           
           // -MULTIPLICITIES:
           //  Secondary multiplicities for all particle types.
           //  Total pion multiplicity, dito separatly for pi+, pi-, pi0
           G4int piMulti = 0;
           G4int piPlusMulti = 0;
           G4int piMinusMulti = 0;
           G4int piZeroMulti = 0;
           //  Total (anti)nucleon multiplicity, dito separatly for p, n
           G4int nucleonMulti = 0;
           G4int protonMulti = 0;
           G4int neutronMulti = 0;
           G4int antiNucleonMulti = 0;
           G4int antiProtonMulti = 0;
           G4int antiNeutronMulti = 0;
           //  Total fragment muliplicity, dito separatly for d, t, alpha
           G4int fragmentMulti = 0;
           G4int deuteronMulti = 0;
           G4int tritonMulti = 0;
           G4int alphaMulti = 0;
           //  Total muliplicity for strange particles,
           //  dito separately for the individual kinds
           G4int kPlusMulti = 0;
           G4int kMinusMulti = 0;
           G4int kZeroLongMulti = 0;
           G4int kZeroShortMulti = 0;
           G4int lambdaMulti = 0;
           G4int antiLambdaMulti = 0;
           G4int sigmaPlusMulti = 0;
           G4int sigmaMinusMulti = 0;
           G4int sigmaZeroMulti = 0;
           G4int xiMinusMulti = 0;
           G4int xiZeroMulti = 0;
           G4int omegaMulti = 0;
           G4int antiSigmaPlusMulti = 0;
           G4int antiSigmaMinusMulti = 0;
           G4int antiSigmaZeroMulti = 0;
           G4int antiXiMinusMulti = 0;
           G4int antiXiZeroMulti = 0;
           G4int antiOmegaMulti = 0;
           G4int strangeMesonMulti = 0;
           G4int strangeBarionMulti = 0;
           G4int strangeAntiBarionMulti = 0;
           //  Total multiplicity of final state particles
           G4int positiveMulti = 0;
           G4int negativeMulti = 0;
           G4int neutralMulti = 0;
           G4int totalMulti = 0;
           // -E-P BALANCE:
           //  Interaction-wise E-P balance
           G4double epBalanceInter = 0;
           // -MOMENTUM DISTRIBUTIONS:
           //  leading final state particle momentum: px, py, pz, E
           //  2nd leading final state particle momentum: px, py, pz, E
           //  3rd leading final state particle momentum: px, py, pz, E
           //  ...
           //  10th leading final state particle momentum: px, py, pz, E
           //  sum of the residual momenta: px, py, pz, E
           G4double sortedPx[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
           G4double restPx = 0;
           G4double sortedPy[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
           G4double restPy = 0;
           G4double sortedPz[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
           G4double restPz = 0;
           G4double sortedEnergy[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
           G4double restEnergy = 0;
           //  Total momentum of pions
           G4double pionPx = 0;
           G4double pionPy = 0;
           G4double pionPz = 0;
           G4double pionEnergy = 0;
           //  Total momentum of Kaons
           G4double kaonPx = 0;
           G4double kaonPy = 0;
           G4double kaonPz = 0;
           G4double kaonEnergy = 0;
           // Total momentum of Nucleons
           G4double nucleonPx = 0;
           G4double nucleonPy = 0;
           G4double nucleonPz = 0;
           G4double nucleonEnergy = 0;
           // fragment energies and momenta
           G4double fragmentPx = 0;
           G4double fragmentPy = 0;
           G4double fragmentPz = 0;
           G4double fragmentEnergy = 0;
           G4double deuteronPx = 0;
           G4double deuteronPy = 0;
           G4double deuteronPz = 0;
           G4double deuteronEnergy = 0;
           G4double tritonPx = 0;
           G4double tritonPy = 0;
           G4double tritonPz = 0;
           G4double tritonEnergy = 0;
           G4double alphaPx = 0;
           G4double alphaPy = 0;
           G4double alphaPz = 0;
           G4double alphaEnergy = 0;
           // Total momentum of pi+, pi-, pi0
           G4double piPlusPx = 0;
           G4double piPlusPy = 0;
           G4double piPlusPz = 0;
           G4double piPlusEnergy = 0;
           G4double piMinusPx = 0;
           G4double piMinusPy = 0;
           G4double piMinusPz = 0;
           G4double piMinusEnergy = 0;
           G4double piZeroPx = 0;
           G4double piZeroPy = 0;
           G4double piZeroPz = 0;
           G4double piZeroEnergy = 0;
           // Total momentum of p, n
           G4double protonPx = 0;
           G4double protonPy = 0;
           G4double protonPz = 0;
           G4double protonEnergy = 0;
           G4double neutronPx = 0;
           G4double neutronPy = 0;
           G4double neutronPz = 0;
           G4double neutronEnergy = 0;
           // Total momentum of K+, K-, Kl, Ks
           G4double kPlusPx = 0;
           G4double kPlusPy = 0;
           G4double kPlusPz = 0;
           G4double kPlusEnergy = 0;
           G4double kMinusPx = 0;
           G4double kMinusPy = 0;
           G4double kMinusPz = 0;
           G4double kMinusEnergy = 0;
           G4double kZeroLongPx = 0;
           G4double kZeroLongPy = 0;
           G4double kZeroLongPz = 0;
           G4double kZeroLongEnergy = 0;
           G4double kZeroShortPx = 0;
           G4double kZeroShortPy = 0;
           G4double kZeroShortPz = 0;
           G4double kZeroShortEnergy = 0;
           // momenta for strange barions
           G4double strangeBarionPx = 0;
           G4double strangeBarionPy = 0;
           G4double strangeBarionPz = 0;
           G4double strangeBarionEnergy = 0;
           G4double strangeAntiBarionPx = 0;
           G4double strangeAntiBarionPy = 0;
           G4double strangeAntiBarionPz = 0;
           G4double strangeAntiBarionEnergy = 0;
           // prepare the analysis current quantities
           G4double currentPx = 0;
           G4double currentPy = 0;
           G4double currentPz = 0;
           G4double currentEnergy= 0;
           G4int currentCharge =0;
           G4ParticleDefinition * currentDefinition;
           
           // now get, and Store the stuff
           
           if( aFinalState->GetStatusChange() == fAlive )
           {
             // include the primary in the analysis
             // later (add it as a secondary)
           }
           for( G4int ii=0; ii<aFinalState->GetNumberOfSecondaries(); ++ii )
           {
             // all the secondaries
             G4Track* currentTrack; 
             currentTrack = aFinalState->GetSecondary(ii);
             currentDefinition = currentTrack->GetDefinition();
             currentEnergy = currentTrack->GetKineticEnergy();
             currentPx = currentTrack->GetMomentum().x();
             currentPy = currentTrack->GetMomentum().y();
             currentPz = currentTrack->GetMomentum().z();
             currentCharge = currentDefinition->GetPDGCharge();
             
             if( currentCharge >  0 )++positiveMulti;
             if( currentCharge <  0 )++negativeMulti;
             if( currentCharge == 0 )++neutralMulti;
             ++totalMulti;
             delete aFinalState->GetSecondary(ii);
             if( currentDefinition == theProton )
             {
               ++nucleonMulti;
               nucleonPx += currentPx;
               nucleonPy += currentPy;
               nucleonPz += currentPz;
               nucleonEnergy += currentEnergy;
               ++protonMulti;
               protonPx += currentPx;
               protonPy += currentPy;
               protonPz += currentPz;
               protonEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  )==0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition==theAntiProton )
             {
               ++antiNucleonMulti;
               ++antiProtonMulti;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  )==0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theNeutron )
             {
               ++nucleonMulti;
               nucleonPx += currentPx;
               nucleonPy += currentPy;
               nucleonPz += currentPz;
               nucleonEnergy += currentEnergy;
               ++neutronMulti;
               neutronPx += currentPx;
               neutronPy += currentPy;
               neutronPz += currentPz;
               neutronEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theAntiNeutron )
             {
               ++antiNucleonMulti;
               ++antiNeutronMulti;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == thePionPlus )
             {
               ++piMulti;
               pionPx += currentPx;
               pionPy += currentPy;
               pionPz += currentPz;
               pionEnergy += currentEnergy;
               ++piPlusMulti;
               piPlusPx += currentPx;
               piPlusPy += currentPy;
               piPlusPz += currentPz;
               piPlusEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == thePionMinus )
             {
               ++piMulti;
               pionPx += currentPx;
               pionPy += currentPy;
               pionPz += currentPz;
               pionEnergy += currentEnergy;
               ++piMinusMulti;
               piMinusPx += currentPx;
               piMinusPy += currentPy;
               piMinusPz += currentPz;
               piMinusEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == G4PionZero::PionZero() )
             {
               ++piMulti;
               pionPx += currentPx;
               pionPy += currentPy;
               pionPz += currentPz;
               pionEnergy += currentEnergy;
               ++piZeroMulti;
               piZeroPx += currentPx;
               piZeroPy += currentPy;
               piZeroPz += currentPz;
               piZeroEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theDeuteron )
             {
               ++fragmentMulti;
               fragmentPx += currentPx;
               fragmentPy += currentPy;
               fragmentPz += currentPz;
               fragmentEnergy += currentEnergy;
               ++deuteronMulti;
               deuteronPx += currentPx;
               deuteronPy += currentPy;
               deuteronPz += currentPz;
               deuteronEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theTriton )
             {
               ++fragmentMulti;
               fragmentPx += currentPx;
               fragmentPy += currentPy;
               fragmentPz += currentPz;
               fragmentEnergy += currentEnergy;
               ++tritonMulti;
               tritonPx += currentPx;
               tritonPy += currentPy;
               tritonPz += currentPz;
               tritonEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theAlpha )
             {
               ++fragmentMulti;
               fragmentPx += currentPx;
               fragmentPy += currentPy;
               fragmentPz += currentPz;
               fragmentEnergy += currentEnergy;
               ++alphaMulti;
               alphaPx += currentPx;
               alphaPy += currentPy;
               alphaPz += currentPz;
               alphaEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theKaonPlus )
             {
               ++strangeMesonMulti;
               ++kPlusMulti;
               kaonPx += currentPx;
               kaonPy += currentPy;
               kaonPz += currentPz;
               kaonEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theKaonMinus )
             {
               ++strangeMesonMulti;
               ++kMinusMulti;
               kaonPx += currentPx;
               kaonPy += currentPy;
               kaonPz += currentPz;
               kaonEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theKaonLong )
             {
               ++strangeMesonMulti;
               ++kZeroLongMulti;
               kaonPx += currentPx;
               kaonPy += currentPy;
               kaonPz += currentPz;
               kaonEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theKaonShort )
             {
               ++strangeMesonMulti;
               ++kZeroShortMulti;
               kaonPx += currentPx;
               kaonPy += currentPy;
               kaonPz += currentPz;
               kaonEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theLambda )
             {
               ++strangeBarionMulti;
               strangeBarionPx += currentPx;
               strangeBarionPy += currentPy;
               strangeBarionPz += currentPz;
               strangeBarionEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theAntiLambda )
             {
               ++strangeAntiBarionMulti;
               strangeAntiBarionPx += currentPx;
               strangeAntiBarionPy += currentPy;
               strangeAntiBarionPz += currentPz;
               strangeAntiBarionEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theSigmaPlus )
             {
               ++strangeBarionMulti;
               strangeBarionPx += currentPx;
               strangeBarionPy += currentPy;
               strangeBarionPz += currentPz;
               strangeBarionEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theAntiSigmaPlus )
             {
               ++strangeAntiBarionMulti;
               strangeAntiBarionPx += currentPx;
               strangeAntiBarionPy += currentPy;
               strangeAntiBarionPz += currentPz;
               strangeAntiBarionEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theSigmaMinus )
             {
               ++strangeBarionMulti;
               strangeBarionPx += currentPx;
               strangeBarionPy += currentPy;
               strangeBarionPz += currentPz;
               strangeBarionEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theAntiSigmaMinus )
             {
               ++strangeAntiBarionMulti;
               strangeAntiBarionPx += currentPx;
               strangeAntiBarionPy += currentPy;
               strangeAntiBarionPz += currentPz;
               strangeAntiBarionEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theXiMinus )
             {
               ++strangeBarionMulti;
               strangeBarionPx += currentPx;
               strangeBarionPy += currentPy;
               strangeBarionPz += currentPz;
               strangeBarionEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theAntiXiMinus )
             {
               ++strangeAntiBarionMulti;
               strangeAntiBarionPx += currentPx;
               strangeAntiBarionPy += currentPy;
               strangeAntiBarionPz += currentPz;
               strangeAntiBarionEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theXiZero )
             {
               ++strangeBarionMulti;
               strangeBarionPx += currentPx;
               strangeBarionPy += currentPy;
               strangeBarionPz += currentPz;
               strangeBarionEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theAntiXiZero )
             {
               ++strangeAntiBarionMulti;
               strangeAntiBarionPx += currentPx;
               strangeAntiBarionPy += currentPy;
               strangeAntiBarionPz += currentPz;
               strangeAntiBarionEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theOmegaMinus )
             {
               ++strangeBarionMulti;
               strangeBarionPx += currentPx;
               strangeBarionPy += currentPy;
               strangeBarionPz += currentPz;
               strangeBarionEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             else if( currentDefinition == theAntiOmegaMinus )
             {
               ++strangeAntiBarionMulti;
               strangeAntiBarionPx += currentPx;
               strangeAntiBarionPy += currentPy;
               strangeAntiBarionPz += currentPz;
               strangeAntiBarionEnergy += currentEnergy;
               epBalanceInter += currentEnergy;
               if( sortEnergies( currentPx, currentPy, currentPz, currentEnergy,
                                 sortedPx,  sortedPy,  sortedPz,  sortedEnergy  ) == 0 )
               {
                 restPx += currentPx;
                 restPy += currentPy;
                 restPz += currentPz;
                 restEnergy += currentEnergy;
               }
             }
             delete currentTrack;
           } // secondary loop
           aFinalState->Clear();
           G4cout << "total multiplicity: "<<totalMulti<<G4endl;
         }  // event loop
       }  // energy loop
     }  // material loop
   }  // particle loop
   return EXIT_SUCCESS;
}


