// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testG4EnergyLossTables.cc,v 1.2 1999-12-15 14:52:04 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//-------------------------------------------------------------------
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"
#include "globals.hh"
#include "G4Timer.hh"
#include "G4MultipleScattering.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4GRSVolume.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4ProcessManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4eEnergyLoss.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4MuEnergyLoss.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4hEnergyLoss.hh"
#include "G4hIonisation.hh"
#include "G4GPILSelection.hh"

#include "G4EnergyLossTables.hh"

// test for G4EnergyLossTables, based on tunemsc
// P. Urban
  
G4VPhysicalVolume* BuildVolume(G4Material* matworld)
//  it builds a simple box filled with material matword .......
{
  G4Box *myWorldBox= new G4Box ("WBox",10000.*cm,10000.*cm,10000.*cm);

  G4LogicalVolume *myWorldLog = new G4LogicalVolume(myWorldBox,matworld,
                                                    "WLog",0,0,0) ;

  G4PVPlacement *myWorldPhys = new G4PVPlacement(0,G4ThreeVector(),
                                                 "WPhys",
                                                 myWorldLog,
                                                 0,false,0) ;

  
  return myWorldPhys ;

}
                                              

int main()
{
  //-------- set output format-------
   G4cout.setf( G4std::ios::scientific, G4std::ios::floatfield );
  //---write results to the file msc.out-----
   G4std::ofstream outFile("msc.out", G4std::ios::out ) ;
   outFile.setf( G4std::ios::scientific, G4std::ios::floatfield );

  //--------- Material definition ---------
  G4Timer theTimer ;
  G4double a, z, ez, density ,temperature,pressure;
  G4State state ;
  G4String name, symbol;
  G4int nel;

  a = 9.012*g/mole;
  density = 1.848*g/cm3;
  G4Material* Be = new G4Material(name="Beryllium", z=4. , a, density);

  a = 12.011*g/mole;
  density = 2.220*g/cm3;
  G4Material* C = new G4Material(name="Carbon", z=6. , a, density);

  a = 26.98*g/mole;
  density = 2.7*g/cm3;
  G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);

  a = 28.09*g/mole;
  density = 2.33*g/cm3;
  G4Material* Si = new G4Material(name="Silicon", z=14., a, density);

  a = 55.85*g/mole;
  density = 7.87*g/cm3;
  G4Material* Fe = new G4Material(name="Iron", z=26., a, density);

  a = 63.540*g/mole;
  density = 8.960*g/cm3;
  G4Material* Cu = new G4Material(name="Copper", z=29., a, density);

  a = 196.97*g/mole;
  density = 19.32*g/cm3;
  G4Material* Au = new G4Material(name="Gold", z=79., a, density);

  a = 207.19*g/mole;
  density = 11.35*g/cm3;
  G4Material* Pb = new G4Material(name="Lead", z=82., a, density);

  a = 238.03*g/mole;
  density = 18.7*g/cm3;
  G4Material* U = new G4Material(name="Uranium", z=92., a, density);

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", ez=7., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxigen", symbol="O", ez=8., a);
  density = 1.29e-03*g/cm3;
  state = kStateGas ;
  temperature = 273.*kelvin ;
  pressure = 1.*atmosphere ;
  G4Material* Air = new G4Material(name="Air", density, nel=2 ,
                                   state ,temperature , pressure ) ;
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);

  a = 0. ;          
  density = 0. ;         
  G4Material* Vac= new G4Material(name="Vacuum",z=0., a, density,kVacuum);

  static const G4MaterialTable* theMaterialTable ;
  G4Material* apttoMaterial ;
  G4String MaterialName ; 

//--------- Particle definition ---------

  G4Gamma* theGamma = G4Gamma::GammaDefinition();

  G4Electron* theElectron = G4Electron::ElectronDefinition();
  G4Positron* thePositron = G4Positron::PositronDefinition();

  G4MuonPlus* theMuonPlus = G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus* theMuonMinus = G4MuonMinus::MuonMinusDefinition();

  G4Proton* theProton = G4Proton::ProtonDefinition();
  G4AntiProton* theAntiProton = G4AntiProton::AntiProtonDefinition();
  G4PionPlus* thePionPlus = G4PionPlus::PionPlusDefinition();
  G4PionMinus* thePionMinus = G4PionMinus::PionMinusDefinition();
  G4KaonPlus* theKaonPlus = G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus* theKaonMinus = G4KaonMinus::KaonMinusDefinition();

//--------- Process definition ---------
  G4eIonisation theElectronIonisation,thePositronIonisation;
  G4eBremsstrahlung theElectronBremsstrahlung,thePositronBremsstrahlung;
//..........e-/e+....................................................
  G4MultipleScattering theElectronMultipleScattering,thePositronMultipleScattering ;
       
  G4ProcessManager* theElectronProcessManager  = theElectron->GetProcessManager();
  G4ProcessManager* thePositronProcessManager  = thePositron->GetProcessManager();
  
  theElectronProcessManager->AddProcess(&theElectronMultipleScattering,-1,0,0) ;
  thePositronProcessManager->AddProcess(&thePositronMultipleScattering,-1,0,0) ;
  theElectronProcessManager->AddProcess(&theElectronIonisation,-1,1,1) ;
  theElectronProcessManager->AddProcess(&theElectronBremsstrahlung,-1,-1,2) ;
  thePositronProcessManager->AddProcess(&thePositronIonisation,-1,1,1) ;
  thePositronProcessManager->AddProcess(&thePositronBremsstrahlung,-1,-1,2) ;

//..........mu+/mu-................................................
  G4MuIonisation theMuonPlusIonisation,theMuonMinusIonisation ;
  G4MuBremsstrahlung theMuonPlusBremsstrahlung,theMuonMinusBremsstrahlung ;
  G4MuPairProduction theMuonPlusPairProduction,theMuonMinusPairProduction ;  
  G4MultipleScattering theMuonPlusMultipleScattering,theMuonMinusMultipleScattering;

  G4ProcessManager* theMuonPlusProcessManager = theMuonPlus->GetProcessManager();
  theMuonPlusProcessManager->AddProcess(&theMuonPlusMultipleScattering,-1,0,0) ;
  theMuonPlusProcessManager->AddProcess(&theMuonPlusIonisation,-1,1,1);
  theMuonPlusProcessManager->AddProcess(&theMuonPlusBremsstrahlung,-1,-1,2);
  theMuonPlusProcessManager->AddProcess(&theMuonPlusPairProduction,-1,-1,3);

  G4ProcessManager* theMuonMinusProcessManager = theMuonMinus->GetProcessManager();
  theMuonMinusProcessManager->AddProcess(&theMuonMinusMultipleScattering,-1,0,0) ;
  theMuonMinusProcessManager->AddProcess(&theMuonMinusIonisation,-1,1,1);
  theMuonMinusProcessManager->AddProcess(&theMuonMinusBremsstrahlung,-1,-1,2);
  theMuonMinusProcessManager->AddProcess(&theMuonMinusPairProduction,-1,-1,3);

//------multiple instantiation of G4MultipleScattering  for hadrons-----------------
  G4MultipleScattering theProtonMultipleScattering,theAntiProtonMultipleScattering;
  G4hIonisation theProtonIonisation,theAntiProtonIonisation ;
  G4MultipleScattering thePionPlusMultipleScattering,thePionMinusMultipleScattering;
  G4hIonisation thePionPlusIonisation,thePionMinusIonisation;
  G4MultipleScattering theKaonPlusMultipleScattering,theKaonMinusMultipleScattering;
  G4hIonisation theKaonPlusIonisation,theKaonMinusIonisation;

  G4ProcessManager* theProtonProcessManager = theProton->GetProcessManager();
  theProtonProcessManager->AddProcess(&theProtonMultipleScattering,-1,0,0) ;
  theProtonProcessManager->AddProcess(&theProtonIonisation,-1,1,1);

  G4ProcessManager* thePionPlusProcessManager = thePionPlus->GetProcessManager();
  thePionPlusProcessManager->AddProcess(&thePionPlusMultipleScattering,-1,0,0) ;
  thePionPlusProcessManager->AddProcess(&thePionPlusIonisation,-1,1,1);

  G4ProcessManager* theKaonPlusProcessManager = theKaonPlus->GetProcessManager();
  theKaonPlusProcessManager->AddProcess(&theKaonPlusMultipleScattering,-1,0,0) ;
  theKaonPlusProcessManager->AddProcess(&theKaonPlusIonisation,-1,1,1);

  G4ProcessManager* theAntiProtonProcessManager = theAntiProton->GetProcessManager();
  theAntiProtonProcessManager->AddProcess(&theAntiProtonMultipleScattering,-1,0,0) ;
  theAntiProtonProcessManager->AddProcess(&theAntiProtonIonisation,-1,1,1);

  G4ProcessManager* thePionMinusProcessManager = thePionMinus->GetProcessManager();
  thePionMinusProcessManager->AddProcess(&thePionMinusMultipleScattering,-1,0,0) ;
  thePionMinusProcessManager->AddProcess(&thePionMinusIonisation,-1,1,1);

  G4ProcessManager* theKaonMinusProcessManager = theKaonMinus->GetProcessManager();
  theKaonMinusProcessManager->AddProcess(&theKaonMinusMultipleScattering,-1,0,0) ;
  theKaonMinusProcessManager->AddProcess(&theKaonMinusIonisation,-1,1,1);

  G4GPILSelection selection;
  G4ParticleWithCuts* theParticle ;
  G4MultipleScattering* theParticleMultipleScattering;
  G4ProcessManager* theParticleProcessManager ;
  G4String confirm ;
  G4int i1e ;
  G4double theta1e,th1,th2 ;
  G4double theta1edata,errth,lambdadata,lambdadatamin,lambdadatamax ;

  // loop until no particle is selected
  while (1) {

  G4cout << "Do you want the electron as particle (yes/no)?" << G4std::flush;
  G4cin >> confirm ;
  if(confirm == "yes")
  {
    theParticle = theElectron ;
    theParticleMultipleScattering=&theElectronMultipleScattering;
    theParticleProcessManager=theElectronProcessManager;
    outFile << " ----------particle = electron -------------" << G4endl;
  }
  else
  {    
    G4cout << "Do you want the positron as particle (yes/no)?" << G4std::flush;
    G4cin >> confirm ;
    if(confirm == "yes")
    {
      theParticle = thePositron ;
      theParticleMultipleScattering=&thePositronMultipleScattering;
      theParticleProcessManager=thePositronProcessManager;
      outFile << " ----------particle = positron -------------" << G4endl;
    }
    else
    {
      G4cout << "Do you want the mu+ as particle (yes/no)?" << G4std::flush;
      G4cin >> confirm ;
      if(confirm == "yes")
      {
        theParticle = theMuonPlus ;
        theParticleMultipleScattering=&theMuonPlusMultipleScattering;
        theParticleProcessManager=theMuonPlusProcessManager;
        outFile << " --------particle = mu+ -------------" << G4endl;
      }
      else
      {
        G4cout << "Do you want the mu- as particle (yes/no)?" << G4std::flush;
        G4cin >> confirm ;
        if(confirm == "yes")
        {
          theParticle = theMuonMinus ;
          theParticleMultipleScattering=&theMuonMinusMultipleScattering;
          theParticleProcessManager=theMuonMinusProcessManager;
          outFile << " --------particle = mu- -------------" << G4endl;
      }
      else
  {
  G4cout << " Do you want the proton as particle (yes/no)? " << G4std::flush;
  G4cin >> confirm ;
  if(confirm == "yes")
  {
    theParticle = theProton;
    theParticleMultipleScattering=&theProtonMultipleScattering;
    theParticleProcessManager=theProtonProcessManager;
    outFile << " ---------- particle = proton ----------------" << G4endl;
  }
  else
  {
     G4cout << " Do you want the antiproton as particle (yes/no)? " << G4std::flush;
     G4cin >> confirm ;
     if(confirm == "yes")
     {
        theParticle = theAntiProton;
        theParticleMultipleScattering=&theAntiProtonMultipleScattering;
        theParticleProcessManager=theAntiProtonProcessManager;
        outFile << " ---------- particle = antiproton ----------------" << G4endl;
     }
     else
     {
      G4cout << " Do you want the pi+ as particle (yes/no)? " << G4std::flush;
      G4cin >> confirm ;
      if(confirm == "yes")
      {
      theParticle = thePionPlus;
      theParticleMultipleScattering=&thePionPlusMultipleScattering;
      theParticleProcessManager=thePionPlusProcessManager;
      outFile << " ---------- particle = pi+ ----------------" << G4endl;
      }
      else
      {
        G4cout << " Do you want the pi- as particle (yes/no)? " << G4std::flush;
        G4cin >> confirm ;
        if(confirm == "yes")
        {
        theParticle = thePionMinus;
        theParticleMultipleScattering=&thePionMinusMultipleScattering;
        theParticleProcessManager=thePionMinusProcessManager;
        outFile << " ---------- particle = pi- ----------------" << G4endl;
        } 
        else
        {
          G4cout << " Do you want the K+ as particle (yes/no)? " << G4std::flush;
          G4cin >> confirm ;
          if(confirm == "yes")
          {
          theParticle = theKaonPlus;
          theParticleMultipleScattering=&theKaonPlusMultipleScattering;
          theParticleProcessManager=theKaonPlusProcessManager;
          outFile << " ---------- particle = K+ ----------------" << G4endl;
          }
          else
          {
            G4cout << " Do you want the K- as particle (yes/no)? " << G4std::flush;
            G4cin >> confirm ;
            if(confirm == "yes")
            {
            theParticle = theKaonMinus;
            theParticleMultipleScattering=&theKaonMinusMultipleScattering;
            theParticleProcessManager=theKaonMinusProcessManager;
            outFile << " ---------- particle = K- ----------------" << G4endl;
            }
            else
            {
             G4cout << " There is no other particle in the test." << G4endl;
             return EXIT_SUCCESS;
            }
          }
        }
       }
      }
     }
    }
   }
  }
 }

  G4ForceCondition cond;
  G4ForceCondition* condition=&cond ;

  G4double* ElectronKineticEnergyCuts ;
  G4double* PositronKineticEnergyCuts ;
  G4double* GammaKineticEnergyCuts ;
  G4double* ParticleKineticEnergyCuts ;

  theMaterialTable = G4Material::GetMaterialTable() ;

  G4double cutinrange ;


  G4cout << "give cuts in range" << G4endl ;

  G4cout << "cut for GAMMA in mm =" ;
  G4cin >> cutinrange ; 
  theGamma->SetCuts(cutinrange) ;

  GammaKineticEnergyCuts = theGamma->GetCutsInEnergy() ;

  G4cout << "cut for ELECTRON in mm =" ;
  G4cin >> cutinrange ; 
  theElectron->SetCuts(cutinrange) ;

  ElectronKineticEnergyCuts = theElectron->GetCutsInEnergy() ;

  G4cout << "cut for POSITRON in mm =" ;
  G4cin >> cutinrange ; 
  thePositron->SetCuts(cutinrange) ;

  PositronKineticEnergyCuts = thePositron->GetCutsInEnergy() ;

//****************************************************
// setcut for the selected particle (if it is not e-/e+)
 if((theParticle != theElectron) && (theParticle != thePositron))
 {
  G4cout << "cut for the selected particle in mm =" ;
  G4cin >> cutinrange ; 
  theParticle->SetCuts(cutinrange) ;

  ParticleKineticEnergyCuts = theParticle->GetEnergyCuts() ;
 }

  G4double energy, momentum, mass;
  G4ProcessVector* palongget ;
  G4ProcessVector* palongdo ;
  G4ProcessVector* ppostget ;
  G4ProcessVector* ppostdo ;

    mass = theParticle->GetPDGMass() ;  
    energy = 1.*GeV + mass ;
    momentum=sqrt(energy*energy-mass*mass) ;
    G4ParticleMomentum theMomentum(momentum,0.,0.);
    G4double pModule = theMomentum.mag();
    G4DynamicParticle aParticle(theParticle,energy,theMomentum);
    aParticle.SetKineticEnergy(energy-mass);
    
    outFile << "  " << G4endl;
    outFile << " M S C test **********************************************" << G4endl ;
    outFile << "  " << G4endl;
    palongget = aParticle.GetDefinition()->GetProcessManager()
                                 ->GetAlongStepProcessVector(0);
    ppostget = aParticle.GetDefinition()->GetProcessManager()
                                 ->GetPostStepProcessVector(0);
    palongdo = aParticle.GetDefinition()->GetProcessManager()
                                 ->GetAlongStepProcessVector(1);
    ppostdo = aParticle.GetDefinition()->GetProcessManager()
                                 ->GetPostStepProcessVector(1);

//---------------------------------- Physics --------------------------------

  G4int itry=1, Ntry=1, Nstart, ir;
  G4double r ;

//**************************************************************************
  const G4int Nbin=100 ;
  G4double TkinMeV[Nbin]  =
              {0.0002, 0.0005, 0.001,0.0015,0.002,0.003,0.004,0.005,0.006,0.008,
               0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.08,
               0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.8,
               1.,1.5,2.,3.,4.,5.,6.,8.,
               10.,15.,20.,30.,40.,50.,60.,80.,
               100.,150.,200.,300.,400.,500.,600.,800.,
               1.0e3,1.5e3,2.0e3,3.0e3,4.0e3,5.0e3,6.0e3,8.0e3,
               1.0e4,1.5e4,2.0e4,3.0e4,4.0e4,5.0e4,6.0e4,8.0e4,
               1.0e5,1.5e5,2.0e5,3.0e5,4.0e5,5.0e5,6.0e5,8.0e5,
               1.0e6,1.5e6,2.0e6,3.0e6,4.0e6,5.0e6,6.0e6,8.0e6,
               1.0e7,1.5e7,2.0e7,3.0e7,4.0e7,5.0e7,6.0e7,8.0e7,
               1.0e8,1.5e8,2.0e8,3.0e8,4.0e8,5.0e8,6.0e8,8.0e8,
               1.0e9, 5.0e9} ; 
         
    G4int J=-1 ;

    G4double lambda,trueStep,geomStep,stepLimit,
             previousStepSize,currentMinimumStep ;
    G4ParticleChange* aParticleChange ;
    
    // material loop
    for (J=0; J<theMaterialTable->length(); J++) {

    apttoMaterial = (*theMaterialTable)[ J ] ;
    MaterialName = apttoMaterial->GetName() ; 
    G4cout << G4endl << "material=" << MaterialName << G4endl << G4endl;

//---------- Volume definition ---------------------

    G4VPhysicalVolume* myVolume ;

    myVolume = BuildVolume(apttoMaterial) ;

//--------- track and Step definition (for this test ONLY!)------------
    G4ThreeVector aPosition(0.,0.,0.);
    const G4ThreeVector aDirection(0.,0.,1.) ;
    const G4ThreeVector transl(0.,0.,0.) ;

    G4double aTime = 0. ;

    G4Track* tracke = new G4Track(&aParticle,aTime,aPosition) ;
    G4Track& trackele = (*tracke) ;
    G4GRSVolume*  touche = new G4GRSVolume(myVolume,NULL,transl);
    (*tracke).SetTouchable(touche);

    (*tracke).SetMomentumDirection(aDirection) ;

    G4Step* Step = new G4Step() ;
  //******************************************
    tracke->SetStep(Step);
    (*Step).InitializeStep(tracke) ;

    const G4Step& Step = (*Step) ;

    G4double CurrentSafety=1000.;
    G4double& currentSafety = CurrentSafety ;

    G4StepPoint* aPoint = new G4StepPoint();
    (*aPoint).SetPosition(aPosition) ;
    G4double safety = 10000.*cm ;
    (*aPoint).SetSafety(safety) ;

    (*Step).SetPostStepPoint(aPoint) ;

    G4int ib,ev ;
    G4double xc,yc,zc,later ; 
//**************************************************************************

 //---------------------------------------------------------------------
    G4double distr[100]=
                       {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    G4double distrx[100]=
                       {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

    G4int nev[100]=
                   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    G4int nevx[100]=
                   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    G4double thetamean,thetamean2,dthetamean,errdistr ;

    outFile << G4endl ;
    outFile << "  " << MaterialName << "  PostStepDoIt (scattering) test " << G4endl ;
    outFile << "  +++++++++++++++++++++++++++++++++++++++++++++++++" << G4endl ;
    outFile << G4endl ;

    G4double wg[100] ;

    G4double over,overx ;
    G4double fwg, costheta,theta ;
    G4int ibin,iw ;
    G4double TMeV ;
    G4double dthetadeg,dtheta ;
    G4int flagdeg ;

//*************************************************************************************
// example: how to use the GetRange() fuctions............

    for (G4int i=0; i<Nbin; i++) {
      G4double TMeV= TkinMeV[i]*MeV;

      G4double rangestart;
      if( theParticle == theElectron )
	rangestart = theElectronIonisation.OldGetRange(theElectron,TMeV,apttoMaterial) ;
      if( theParticle == thePositron )
	rangestart = thePositronIonisation.OldGetRange(thePositron,TMeV,apttoMaterial) ;
      if( theParticle == theMuonPlus )
	rangestart = theMuonPlusIonisation.OldGetRange(theMuonPlus,TMeV,apttoMaterial) ;
      if( theParticle == theMuonMinus )
	rangestart = theMuonMinusIonisation.OldGetRange(theMuonMinus,TMeV,apttoMaterial) ; 
      if( theParticle == theProton )
	rangestart = theProtonIonisation.OldGetRange(theProton,TMeV,apttoMaterial) ;
      if( theParticle == theAntiProton )
	rangestart = theAntiProtonIonisation.OldGetRange(theAntiProton,TMeV,apttoMaterial) ; 
      if( theParticle == thePionPlus )
	rangestart = thePionPlusIonisation.OldGetRange(thePionPlus,TMeV,apttoMaterial) ; 
      if( theParticle == thePionMinus )
	rangestart = thePionMinusIonisation.OldGetRange(thePionMinus,TMeV,apttoMaterial) ; 
      if( theParticle == theKaonPlus )
	rangestart = theKaonPlusIonisation.OldGetRange(theKaonPlus,TMeV,apttoMaterial) ; 
      if( theParticle == theKaonMinus )
	rangestart = theKaonMinusIonisation.OldGetRange(theKaonMinus,TMeV,apttoMaterial) ; 
      G4cout << " energy (MeV) = " << TMeV << G4endl;
      G4cout << " range (OldGetRange) = " << rangestart << G4endl ;   

      G4double range2= G4EnergyLossTables::GetRange(theParticle, TMeV, apttoMaterial);

      G4cout << " range (G4EnergyLossTables::GetRange) = " << range2 << G4endl ;   
    }
            
    } // of material loop

  } // of largest loop

  return EXIT_SUCCESS;
}
