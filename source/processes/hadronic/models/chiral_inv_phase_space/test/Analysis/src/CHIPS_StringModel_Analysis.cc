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
// $Id: CHIPS_StringModel_Analysis.cc,v 1.14 2009-12-15 10:13:35 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Created by Johannes Peter Wellisch, 22 Apr 1997: full test-suite coded
// Updated and developed by Mikhail Kosov, 9 Dec 2009: for physics lists
    
using namespace std;

#include "ParticleInfo.h"

#include "G4ios.hh"
#include <fstream>
#include <iomanip>

#include "G4Timer.hh"
 
#include "G4Material.hh"
#include "G4Navigator.hh"
#include "G4RunManager.hh"
#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "G4UserRunAction.hh"
#include "G4UserTrackingAction.hh"

//---LOOKHERE---
// Uncomment these includes files only if you use a Geant4 version before 9.2 .
#include "LHEP.hh"
#include "QGSP_BERT.hh"
#include "QGSC_BERT.hh"
#include "FTF_BIC.hh"
#include "CHIPS.hh"


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
#include "G4StringChipsParticleLevelInterface.hh"
#include "G4StringChipsInterface.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4StringModel.hh"
#include "G4VPreCompoundModel.hh"
#include "G4FTFModel.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4LundStringFragmentation.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4QGSParticipants.hh"

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

class TrackingAction : public G4UserTrackingAction
{
public:

  TrackingAction(){;}

  ~TrackingAction(){;}

  void PreUserTrackingAction( const G4Track* ){;}

  void PostUserTrackingAction( const G4Track* ){;}

};

class RunAction: public G4UserRunAction
{
public:

  RunAction(){;}
  virtual ~RunAction(){;}
      
  virtual void BeginOfRunAction( const G4Run* ){;}    
  virtual void EndOfRunAction( const G4Run* ){;}    

};

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:

  PrimaryGeneratorAction(){;}
  ~PrimaryGeneratorAction(){;}

public:

  void GeneratePrimaries( G4Event* ){;}

private:

  G4ParticleGun* particleGun;
};

class DetectorConstruction : public G4VUserDetectorConstruction
{

public:

  DetectorConstruction(){;}
  ~DetectorConstruction(){;}
  
  G4VPhysicalVolume* Construct()
  {
    G4Material* theCu = new G4Material("Copper1", 8.96*g/cm3, 1);
    G4Element* elCu = new G4Element("Copper1", "Cu1", 29., 63.55*g/mole);
    theCu->AddElement( elCu, 1 );
    G4Box* sv = new G4Box ( "SV1",1000*m, 1000*m, 1000*m );
    G4LogicalVolume* lv = new G4LogicalVolume(sv, theCu, "LV1", 0, 0, 0 );
    G4PVPlacement* pv = new G4PVPlacement(0, G4ThreeVector(), "PV1", lv, 0, false, 0 );
    return pv;
  }

  void SetMagField( const G4double){;}
  void SetAbsorberMaterial( const G4String ){;}
  void SetActiveMaterial( const G4String){;}
  // Use by the messenger.

  //G4Material* GetAbsorberMaterial() const{ return new G4Material;}
  //G4Material* GetActiveMaterial() const{ return new G4Material;}

  void SetIsCalHomogeneous( const G4bool ){;}
  void SetIsUnitInLambda( const G4bool ){;}
  void SetAbsorberTotalLength( const G4double ){;}
  void SetCalorimeterRadius( const G4double ){;}
  void SetActiveLayerNumber( const G4int ){;}
  void SetActiveLayerSize( const G4double ){;}
  void SetReadoutLayerNumber( const G4int ){;}
  // To define the calorimeter geometry.

  void SetIsRadiusUnitInLambda( const G4bool ){;}
  void SetRadiusBinSize( const G4double ){;}
  void SetRadiusBinNumber( const G4int ){;}
  // To define the transverse shower analysis.
  
  void UpdateGeometry(){;}

private:

  void DefineMaterials(){;}
  // Define all the materials.

  G4VPhysicalVolume* ConstructCalorimeter()
  {
    G4Material* theCu = new G4Material("Copper2", 8.96*g/cm3, 1);
    G4Element* elCu = new G4Element("Copper2", "Cu2", 29., 63.55*g/mole);
    theCu->AddElement( elCu, 1 );
    G4Box* sv = new G4Box ( "SV2",100*m, 100*m, 100*m );
    G4LogicalVolume* lv = new G4LogicalVolume(sv, theCu, "LV2", 0, 0, 0 );
    G4PVPlacement* pv = new G4PVPlacement(0, G4ThreeVector(), "PV2", lv, 0, false, 0 );
    return pv;
  }
  // To be invoked each time the geometry needs to be updated.

  G4bool areParametersOK(){ return true;}
  // Return true if all the parameters are sensible, false otherwise.

  void PrintParameters(){;}
  // Print the various parameters which define the calorimeter.
};

class SteppingAction : public G4UserSteppingAction
{
public:

  SteppingAction(){;}
  virtual ~SteppingAction(){;}

public:

  virtual void UserSteppingAction( const G4Step* ){;}
  // The main method to define.

  G4double getTotalEdepAllParticles() const {return 0.;}
  // Total deposited energy, in the same event, due to all particles,
  // inside the calorimeter.

  G4int getPrimaryParticleId() const {return 0;}
  G4double getPrimaryParticleEnergy() const {return 0.;}
  // Id and Energy of the primary particle.

  void reset(){;}
  // This method should be invoked at the end of the event,
  // to reset to zero the total deposit energy variables.

private:
 
  //G4double totalEdepAllParticles;
  //G4int    primaryParticleId;
  //G4double primaryParticleEnergy;
  //G4bool   isFirstStepOfTheEvent;
};
class EventAction: public G4UserEventAction
{
public:
  EventAction(){;}
  ~EventAction(){;}

  virtual void BeginOfEventAction( const G4Event* ){;}    
  virtual void EndOfEventAction( const G4Event* ){;}    

private:

  void instanciateSteppingAction(){;}      
  SteppingAction* theSteppingAction;
  //G4Timer* eventTimer;
};

class PiToKPoint
{
 public:
    
   PiToKPoint(double low, double high)
   {
     thePip = thePim = theKp = theKm = 0;
     highPl = high;
     lowPl = low;
   }
   
   bool push_back(const G4DynamicParticle * aPart)
   {
     double pl = aPart->GetMomentum().z();
     if( !(pl <= highPl && pl > lowPl) ) return false;
     
     const G4ParticleDefinition * currentDef = aPart->GetDefinition();
     if     ( currentDef == G4PionPlus::PionPlusDefinition()   ) thePip++;
     else if( currentDef == G4PionMinus::PionMinusDefinition() ) thePim++;
     else if( currentDef == G4KaonPlus::KaonPlusDefinition()   ) theKp++;
     else if( currentDef == G4KaonMinus::KaonMinusDefinition() ) theKm++;
     else  return false;
     return true;
   }
   
   void dump()
   {
     cerr<<"K/PiDump: "<<(highPl+lowPl)/2<<" "<<thePip<<" "<<thePim<<" "<<theKp<<" "<<theKm
         <<" "<<double(theKp)/max(1,thePip)<<" "<<double(theKm)/max(1,thePim)<<endl;
   }
      
 private:   
   double highPl;
   double lowPl;
   int thePip;
   int thePim;
   int theKp;
   int theKm;
   
};

class KtoPi
{
 public:  
   KtoPi(double highCosTh)
   {
     theHighCosTh = highCosTh;
     
     PiToKPoint* p1 = new PiToKPoint(6.8*GeV, 7.2*GeV);
     PiToKPoint* p2 = new PiToKPoint(9.8*GeV, 10.2*GeV);
     PiToKPoint* p3 = new PiToKPoint(14.8*GeV, 15.2*GeV);
     PiToKPoint* p4 = new PiToKPoint(19.6*GeV, 20.4*GeV);
     PiToKPoint* p5 = new PiToKPoint(29.4*GeV, 30.6*GeV);
     PiToKPoint* p6 = new PiToKPoint(39.2*GeV, 40.8*GeV);
     PiToKPoint* p7 = new PiToKPoint(64.*GeV, 66.*GeV);
     PiToKPoint* p8 = new PiToKPoint(132.*GeV, 138.*GeV);
     
     theData.push_back(p1);
     theData.push_back(p2);
     theData.push_back(p3);
     theData.push_back(p4);
     theData.push_back(p5);
     theData.push_back(p6);
     theData.push_back(p7);
     theData.push_back(p8);
   }
   
   void push_back(const G4DynamicParticle * aPart)
   {
     double pl = aPart->GetMomentum().z();
     double pp = aPart->GetMomentum().mag();
     if(pl<6.*GeV) return;
     if(pl/pp<theHighCosTh) return;
     
     for(unsigned i=0; i<theData.size(); ++i) theData[i]->push_back(aPart);
   }
   
   void dump() {for(unsigned i=0; i<theData.size(); ++i) theData[i]->dump();}
   
 private:
   KtoPi(){}
   KtoPi(const KtoPi &){}
 
 private:
   vector<PiToKPoint *> theData;
   double theHighCosTh;
};

int main()
{    
  //RanecuEngine theEngine;
  //HepRandom::setTheEngine( &theEngine );
  //theEngine.showStatus();
  //cout << endl;

  //------ Physics list definition -------
  G4RunManager* runManager = new G4RunManager;
  runManager->SetUserInitialization( new DetectorConstruction ); 
  G4VModularPhysicsList* physList = 0;
  G4int choice=0;
  do
  {
    cout<<"Please select PhysList (1-CHIPS, 2-QGSP, 3-QGSC, 4-FTFB, 5-LHEP): "<<std::flush;
    cin >> choice;
    if      (choice==1) physList = new CHIPS();
    else if (choice==2) physList = new QGSP_BERT();  // @@ DELETE !!
    else if (choice==3) physList = new QGSC_BERT();
    else if (choice==4) physList = new FTF_BIC();
    else if (choice==5) physList = new LHEP();
  } while(choice < 0  || choice > 5);
  runManager->SetUserInitialization( physList );

  // ----- Materials/Elements/Isotopes definition -----
  G4String name, symbol;
  G4double a, iz, density;
  G4int nEl;
    
  G4Material* theCu = new G4Material(name="Copper", density=8.96*g/cm3, nEl=1);
  G4Element* elCu = new G4Element(name="Copper", symbol="Cu", iz=29., a=63.55*g/mole);
  theCu->AddElement( elCu, 1 );
  G4Material* thePb = new G4Material(name="Lead", density=11.35*g/cm3, nEl=1);
  G4Element* elPb = new G4Element(name="Lead", symbol="Pb", iz=82., a=207.19*g/mole);
  thePb->AddElement( elPb, 1 );
  G4Material* theFe = new G4Material(name="Iron", density=7.87*g/cm3, nEl=1);
  G4Element* elFe = new G4Element(name="Iron", symbol="Fe", iz=26., a=55.85*g/mole);
  theFe->AddElement( elFe, 1 );
  G4Material* theW  = new G4Material(name="Tungsten", density=19.30*g/cm3, nEl=1);
  G4Element* elW  = new G4Element(name="Tungston", symbol="W", iz=74., a=183.85*g/mole);
  theW->AddElement( elW, 1 );
  G4Material* theLAr= new G4Material(name="LArgon", density=1.393*g/cm3, nEl=1);
  G4Element* elAr  = new G4Element(name="Argon", symbol="Ar", iz=18., a=39.95*g/mole);
  theLAr->AddElement( elAr, 1 );
  G4Material* thePS = new G4Material(name="PolyStyrene", density=1.032*g/cm3, nEl=2);
  G4Element* elC = new G4Element(name="Carbon", symbol="C", iz=6., a=12.01*g/mole);
  G4Element* elH = new G4Element(name="Hydrogen", symbol="H", iz=1., a=1.01*g/mole);
  thePS->AddElement( elC, 8 );
  thePS->AddElement( elH, 8 );
  G4Material* thePbWO4 = new G4Material(name="PbWO4", density=12.0*g/cm3, nEl=3);
  // approximative number
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", iz=8., a=15.9994*g/mole);
  thePbWO4->AddElement( elPb, 1 );
  thePbWO4->AddElement( elW,  1 );
  thePbWO4->AddElement( elO,  4 );
  // approximate numbers for O
  G4Material* theO = new G4Material(name="Oxygen", density=1.1*g/cm3, nEl=1);
  theO->AddElement( elO,  1 );
  G4Material* theBe = new G4Material(name="Beryllium", density=1.848*g/cm3, nEl=1);
  G4Element* elBe  = new G4Element(name="Beryllium", symbol="Be", iz=4., a=9.01*g/mole);
  theBe->AddElement( elBe, 1 );
  G4Material* theAl = new G4Material(name="Aluminium", density=2.70*g/cm3, nEl=1);
  G4Element* elAl  = new G4Element(name="Aluminium", symbol="Al", iz=13., a=26.98*g/mole);
  theAl->AddElement( elAl, 1 );
  G4Material* theU = new G4Material(name="Uranium", density=18.95*g/cm3, nEl=1);
  G4Element* elU  = new G4Element(name="Uranium", symbol="U", iz=92., a=238.03*g/mole);
  theU->AddElement( elU, 1 );
  G4Material* theBGO = new G4Material(name="BGO", density=2.15*g/cm3, nEl=3);
  G4Element* elBi = new G4Element(name="Bismuth", symbol="Bi", iz=83., a=208.98*g/mole);
  G4Element* elGe = new G4Element(name="Germanium", symbol="Ge", iz=32., a=72.59*g/mole);
  theBGO->AddElement( elBi, 4 );
  theBGO->AddElement( elGe, 3 );
  theBGO->AddElement( elO, 12 );
  G4Material* theNaI = new G4Material(name="NaI", density=3.67*g/cm3, nEl=2);
  G4Element* elNa = new G4Element(name="Sodium", symbol="Na", iz=11., a=22.990*g/mole);
  G4Element* elI = new G4Element(name="Iodine", symbol="I", iz=53., a=126.904*g/mole);
  theNaI->AddElement( elNa, 1 );
  theNaI->AddElement( elI, 1 );
  G4Material *theCsI = new G4Material(name="CsI", density=4.53*g/cm3, nEl=2);
  G4Element* elCs = new G4Element(name="Cesium", symbol="Cs", iz=55., a=132.905*g/mole);
  theCsI->AddElement( elCs, 1 );
  theCsI->AddElement( elI, 1 );
  G4Material* theKapton = new G4Material(name="Kapton", density=1.53*g/cm3, nEl=4); 
  // formula: private communications, see mail.
  theKapton->AddElement( elC, 22 );
  theKapton->AddElement( elH, 10 );
  theKapton->AddElement( elO,  5 );
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a=14.007*g/mole);
  theKapton->AddElement( elN, 2 );
  G4Material* theH = new G4Material(name="Hydrogen", density=1.53*g/cm3, nEl=1); 
  theH->AddElement( elH, 1 );
  G4Material* theC = new G4Material(name="Carbon", density=1.032*g/cm3, nEl=1);
  theC->AddElement( elC, 1 );
  G4Material* theLi6 = new G4Material(name="Li6", density=2.70*g/cm3, nEl=1);
  G4int nIso;
  G4Element* elLi6  = new G4Element(name="Li6", symbol="Li", nIso = 1);
  G4Isotope* isoLi6 = new G4Isotope(name="Li6", 3, 6, a=6.*g/mole);
  elLi6->AddIsotope(isoLi6, 1);
  theLi6->AddElement(elLi6 , 1 );
  G4Material* theTa = new G4Material(name="Tantalum", density=1.53*g/cm3, nEl=1); 
  G4Element* elTa = new G4Element(name="Tantalum", symbol="Ta", iz=73., a=181*g/mole);
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
    
  G4LogicalVolume* LogicalFrame = new G4LogicalVolume(theFrame, (*theMaterialTable)[imat],
                                                      "LFrame", 0, 0, 0 );
    
  G4PVPlacement* PhysicalFrame = new G4PVPlacement(0, G4ThreeVector(), "PFrame",
                                                   LogicalFrame, 0, false, 0 );
  G4RotationMatrix theNull;
  G4ThreeVector theCenter(0,0,0);
  //G4GRSVolume* theTouchable = new G4GRSVolume(PhysicalFrame, &theNull, theCenter);
  G4Navigator* nav = new G4Navigator;
  nav->SetWorldVolume(PhysicalFrame);
  G4TouchableHandle touch(nav->CreateTouchableHistory());

  // ----------- now get all particles of interest --------- 
  const G4int numberOfParticles = 25;
  G4ParticleDefinition* theParticles[numberOfParticles];
  G4ParticleDefinition* theProton = G4Proton::ProtonDefinition();
  cout<<"ProtonProcessManagerPointer="<<theProton->GetProcessManager()<<endl;
  theParticles[ 0] = theProton;
  G4ParticleDefinition* theAntiProton = G4AntiProton::AntiProtonDefinition();
  theParticles[ 1] = theAntiProton;
  G4ParticleDefinition* theNeutron = G4Neutron::NeutronDefinition();
  theParticles[ 2] = theNeutron;
  G4ParticleDefinition* theAntiNeutron = G4AntiNeutron::AntiNeutronDefinition();
  theParticles[ 3] = theAntiNeutron;
  G4ParticleDefinition* thePionPlus = G4PionPlus::PionPlusDefinition();
  theParticles[ 4] = thePionPlus;
  G4ParticleDefinition* thePionMinus = G4PionMinus::PionMinusDefinition();
  theParticles[ 5] = thePionMinus;
  G4ParticleDefinition* theDeuteron = G4Deuteron::DeuteronDefinition();
  theParticles[ 6] = theDeuteron;
  G4ParticleDefinition* theTriton = G4Triton::TritonDefinition();
  theParticles[ 7] = theTriton;
  G4ParticleDefinition* theAlpha = G4Alpha::AlphaDefinition();
  theParticles[ 8] = theAlpha;
  G4ParticleDefinition* theKaonPlus = G4KaonPlus::KaonPlusDefinition();
  theParticles[ 9] = theKaonPlus;
  G4ParticleDefinition* theKaonMinus = G4KaonMinus::KaonMinusDefinition();
  theParticles[10] = theKaonMinus;
  G4ParticleDefinition* theKaonLong = G4KaonZeroLong::KaonZeroLongDefinition();
  theParticles[11] = theKaonLong;
  G4ParticleDefinition* theKaonShort = G4KaonZeroShort::KaonZeroShortDefinition();
  theParticles[12] = theKaonShort;
  G4ParticleDefinition* theLambda = G4Lambda::LambdaDefinition();
  theParticles[13] = theLambda;
  G4ParticleDefinition* theAntiLambda = G4AntiLambda::AntiLambdaDefinition();
  theParticles[14] = theAntiLambda;
  G4ParticleDefinition* theSigmaPlus = G4SigmaPlus::SigmaPlusDefinition();
  theParticles[15] = theSigmaPlus;
  G4ParticleDefinition* theAntiSigmaPlus = G4AntiSigmaPlus::AntiSigmaPlusDefinition();
  theParticles[16] = theAntiSigmaPlus;
  G4ParticleDefinition* theSigmaMinus = G4SigmaMinus::SigmaMinusDefinition();
  theParticles[17] = theSigmaMinus;
  G4ParticleDefinition* theAntiSigmaMinus = G4AntiSigmaMinus::AntiSigmaMinusDefinition();
  theParticles[18] = theAntiSigmaMinus;
  G4ParticleDefinition* theXiMinus = G4XiMinus::XiMinusDefinition();
  theParticles[19] = theXiMinus;
  G4ParticleDefinition* theAntiXiMinus = G4AntiXiMinus::AntiXiMinusDefinition();
  theParticles[20] = theAntiXiMinus;
  G4ParticleDefinition* theXiZero = G4XiZero::XiZeroDefinition();
  theParticles[21] = theXiZero;
  G4ParticleDefinition* theAntiXiZero = G4AntiXiZero::AntiXiZeroDefinition();
  theParticles[22] = theAntiXiZero;
  G4ParticleDefinition* theOmegaMinus = G4OmegaMinus::OmegaMinusDefinition();
  theParticles[23] = theOmegaMinus;
  G4ParticleDefinition* theAntiOmegaMinus = G4AntiOmegaMinus::AntiOmegaMinusDefinition();
  theParticles[24] = theAntiOmegaMinus;
  //------ all the particles are Done ----------    
  
  // Initialize the runManager and define the particle of interest   
  runManager->SetUserAction( new PrimaryGeneratorAction ); 
  runManager->SetUserAction( new RunAction ); 
  runManager->SetUserAction( new EventAction ); 
  runManager->SetUserAction( new TrackingAction ); 
  runManager->Initialize();
  G4int ir=-1; // ir is an index of the projectile particle under the test
  do
  {
    cout << " 0) proton       1) anti_proton  2) neutron      3) anti_neutron" << endl;
    cout << " 4) pi+          5) pi-          6) deuteron     7) triton" << endl;
    cout << " 8) alpha        9) kaon+       10) kaon-       11) kaon0L" << endl;
    cout << "12) kaon0S      13) lambda      14) anti_lambda 15) sigma+" << endl;
    cout << "16) anti_sigma+ 17) sigma-      18) anti_sigma- 19) xi-" << endl;
    cout << "20) anti_xi-    21) xi0         22) anti_xi0    23) omega-" << endl;
    cout << "24) anti_omega-" << endl;
    cout << "Please enter the projectile particle code: "<< std::flush;
    cin >> ir;
  } while( ir < 0 || ir >= numberOfParticles );
  G4ParticleDefinition* theProjectile=theParticles[ir];

  if(theProjectile == G4Proton::ProtonDefinition()) cout<<"Yes. Proton."<<endl;
  G4ProcessManager* theProcessManager = theProjectile->GetProcessManager();
  cout<<"Proj: "<<theProjectile->GetParticleName()<<", ProcMan="<<theProcessManager<<endl;
  G4ProcessVector* theProcessVector = theProcessManager->GetPostStepProcessVector();
  G4int nProc=theProcessVector->size();
  cout<<"For "<<theProjectile->GetParticleName()<<" "<<nProc<<" processes found:"<<endl;
  for(int ip=0; ip<nProc; ++ip)
    cout<<ip<<". "<<(*theProcessVector)[ip]->GetProcessName()<<endl;
  G4int prind=0;
  do
  {
    cout<<"Please select the process index: "<<std::flush;
    cin >> prind;
  } while(prind < 0  || prind >= nProc);
  G4VProcess* theProcess = (*theProcessVector)[prind];
 
  //return EXIT_SUCCESS; // Temporary for debugging
  
  G4ForceCondition* condition = new G4ForceCondition;
  *condition = NotForced;
  
  G4ParticleMomentum theDirection( 0., 0., 1. );
  G4ThreeVector aPosition( 0., 0., 0. );
  G4double aTime = 0.0;
  
  G4StepPoint aStepPoint;
  G4Step aStep;
  aStep.SetPreStepPoint(&aStepPoint);
  //G4double meanFreePath;
  G4double incomingEnergy;
  //G4ParticleDefinition* currentDefinition;
  //G4Track* currentTrack;
  
  G4int kr=-1; // kr = material number, 
  do
  {
    G4cout << " 0) Copper      1) Lead         2) Iron      3) Tungsten" << G4endl;
    G4cout << " 4) LArgon      5) PolyStyrene  6) PbWO4     7) Oxygen" << G4endl;
    G4cout << " 8) Beryllium   9) Aluminium   10) Uranium  11) BGO" << G4endl;
    G4cout << "12) NaI        13) CsI         14) Kapton   15) Hydrogen" << G4endl;
    G4cout << "16) Carbon     17) 6_3_Lithium 18) Tantalum "<< G4endl;
    G4cout << "Please enter the material code: "<<std::flush;
    G4cin >> kr;
  } while( kr < 0 || kr >= numberOfMaterials);
  G4Material* theMaterial=theMaterials[kr];

  G4cout << "Please enter the initial kinetic energy (GeV): "<< std::flush;
  G4cin >> incomingEnergy;
  incomingEnergy *= GeV/MeV;
    
  G4int nEvents;
  G4cout << "Please enter the number of events: "<< std::flush;
  G4cin >> nEvents;
    
  G4int nCPU = 0;
  G4double dCPU = 0.0;
  G4Timer timerEvent;
  G4Timer timerTotal;
  timerTotal.Start();
  
  // Prepare the analysis
  G4String file("../logs/liste.");
  G4String fileName;  // e.g. p_al
  G4double weight=0.; // sigme/A
  if     (kr == 17) 
  {
    fileName = "p_li6";
    weight = 128.*millibarn/6.;
  }
  else if(kr == 8) 
  {
    fileName = "p_be";
    weight = 204.*millibarn/9.;
  }
  else if(kr == 16) 
  {
    fileName = "p_c";
    weight = 243.*millibarn/12.;
  }
  else if(kr == 9) 
  {
    fileName = "p_al";
    weight = 455.*millibarn/27.;
  }
  else if(kr == 0) 
  {
    fileName = "p_cu";
    weight = 823.*millibarn/65.;
  }
  else if(kr == 18) 
  {
    fileName = "p_ta";
    weight = 1690.*millibarn/181.;
  }
  else cout<<"***!*** unknown interaction XS weight: kr="<<kr<<endl;
  G4String modelName; // e.g. _chips
  if     (choice==1) modelName="_chips";
  else if(choice==2) modelName="_qgspb";
  else if(choice==3) modelName="_qgscb";
  else if(choice==4) modelName="_ftfb";
  else if(choice==5) modelName="_lhep";
  else cout<<"***!*** unknown physics list: "<<choice<<endl;
  G4cout<<"===> Running "<<fileName+modelName<<G4endl;
  G4String it = file+fileName+modelName;
  ANAParticleInfo theInformation(weight, it);
  KtoPi theKtoPi(.9);
  // some extra stuff for pi/K ratios @@@
 
  LogicalFrame->SetMaterial( theMaterial ); 
  // --------- Test the PostStepDoIt now  --------------
  G4ParticleChange *aFinalState;
  LogicalFrame->SetMaterial( theMaterial ); 
  G4DynamicParticle aParticle( theProjectile, theDirection, incomingEnergy );
  for(int l=0; l<nEvents; ++l )
  {
    aParticle.SetDefinition( theProjectile );
    aParticle.SetMomentumDirection( theDirection );
    aParticle.SetKineticEnergy( incomingEnergy );
    
    G4Track *aTrack = new G4Track( &aParticle, aTime, aPosition );
    aTrack->SetTouchableHandle(touch); // Set Box touchable history
    aTrack->SetStep( &aStep );
    aStep.SetTrack( aTrack );
    aStepPoint.SetMaterial( theMaterial );
    aStep.SetPreStepPoint(&aStepPoint);
    aStep.SetPostStepPoint(&aStepPoint);
    //cout<<"Event#"<<std::setw(4)<<l<<", Material: "<<theMaterial->GetName()
    //    <<", Particle: "<<theProjectile->GetParticleName()<< endl;
    timerEvent.Start();
    theProcess->PostStepGetPhysicalInteractionLength(*aTrack, 0.1, condition);
    aFinalState = (G4ParticleChange*) (theProcess->PostStepDoIt(*aTrack, aStep) );
    timerEvent.Stop();

    G4int nSec = aFinalState->GetNumberOfSecondaries();
    
    // prepare the analysis current quantities
    
    //cout<<"NUMBER OF SECONDARIES="<<aFinalState->GetNumberOfSecondaries()<<G4endl;
    G4Track* second;
    G4DynamicParticle* aSec;
    for(int isec=0; isec < nSec; ++isec)
    {
      second = aFinalState->GetSecondary(isec);
      aSec = const_cast<G4DynamicParticle *>(second->GetDynamicParticle());
      theKtoPi.push_back(aSec);
      double chg=aSec->GetDefinition()->GetPDGCharge();
      double bar=aSec->GetDefinition()->GetBaryonNumber();
      //cout<<"SECONDAR%# "<<isec<<", charge="<<chg<<", BN="<<bar<<", PDGCode=";
      //if(aSec->GetDefinition()->GetPDGEncoding())
      //     cout<<aSec->GetDefinition()->GetPDGEncoding();
      //else cout<<(int)(1000000000+chg*1000+bar);
      //cout <<", E="<<aSec->GetTotalEnergy()<<", P= "<<aSec->GetMomentum()<<endl;
      // analysis part;
      int theCHIPSCode = aSec->GetDefinition()->GetPDGEncoding();
      if(theCHIPSCode == 0) theCHIPSCode = (G4int)(1000000000+chg*1000+bar);
      ANAParticle aPart(theCHIPSCode, aSec->GetMomentum().x(), aSec->GetMomentum().y(),
                        aSec->GetMomentum().z(), aSec->GetTotalEnergy());
      theInformation.ProcessOne(aPart);
      delete aSec;
      delete second;
    }
    delete aTrack;
    aFinalState->Clear();
    if(l && !(l%1000)) 
    {
      theInformation.Analyse();
      theInformation.Plot(fileName, l);
      theKtoPi.dump();
    }
  } // end of the event loop
  timerTotal.Stop();    
  cout<<"Terminating successfully"<<endl;
  cout<<"TotalTime="<<timerTotal<<", PureTime="<<dCPU<<", #ofEvents:"<<nCPU<<endl; 
  delete condition;
  delete runManager; 
  return EXIT_SUCCESS;
}
/* end of code */
