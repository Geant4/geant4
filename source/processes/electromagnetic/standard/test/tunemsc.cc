// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: tunemsc.cc,v 1.1 1999-01-08 16:32:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//-------------------------------------------------------------------
#include "G4ios.hh"
#include <fstream.h>
#include <iomanip.h>
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

//   It tests the G4MultipleScattering process -------------------------
//   (i.e. the e+/e- , mu+/mu- , h+/h- multiple scattering processes)
//    created by L.Urban on 15/10/97 ----------------------------------
  
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
   G4cout.setf( ios::scientific, ios::floatfield );
  //---write results to the file msc.out-----
   ofstream outFile("msc.out", ios::out ) ;
   outFile.setf( ios::scientific, ios::floatfield );

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

  NEWPARTICLE: ;

  G4cout << "Do you want the electron as particle (yes/no)?" << flush;
  cin >> confirm ;
  if(confirm == "yes")
  {
    theParticle = theElectron ;
    theParticleMultipleScattering=&theElectronMultipleScattering;
    theParticleProcessManager=theElectronProcessManager;
    outFile << " ----------particle = electron -------------" << endl;
  }
  else
  {    
    G4cout << "Do you want the positron as particle (yes/no)?" << flush;
    cin >> confirm ;
    if(confirm == "yes")
    {
      theParticle = thePositron ;
      theParticleMultipleScattering=&thePositronMultipleScattering;
      theParticleProcessManager=thePositronProcessManager;
      outFile << " ----------particle = positron -------------" << endl;
    }
    else
    {
      G4cout << "Do you want the mu+ as particle (yes/no)?" << flush;
      cin >> confirm ;
      if(confirm == "yes")
      {
        theParticle = theMuonPlus ;
        theParticleMultipleScattering=&theMuonPlusMultipleScattering;
        theParticleProcessManager=theMuonPlusProcessManager;
        outFile << " --------particle = mu+ -------------" << endl;
      }
      else
      {
        G4cout << "Do you want the mu- as particle (yes/no)?" << flush;
        cin >> confirm ;
        if(confirm == "yes")
        {
          theParticle = theMuonMinus ;
          theParticleMultipleScattering=&theMuonMinusMultipleScattering;
          theParticleProcessManager=theMuonMinusProcessManager;
          outFile << " --------particle = mu- -------------" << endl;
      }
      else
  {
  G4cout << " Do you want the proton as particle (yes/no)? " << flush;
  cin >> confirm ;
  if(confirm == "yes")
  {
    theParticle = theProton;
    theParticleMultipleScattering=&theProtonMultipleScattering;
    theParticleProcessManager=theProtonProcessManager;
    outFile << " ---------- particle = proton ----------------" << endl;
  }
  else
  {
     G4cout << " Do you want the antiproton as particle (yes/no)? " << flush;
     cin >> confirm ;
     if(confirm == "yes")
     {
        theParticle = theAntiProton;
        theParticleMultipleScattering=&theAntiProtonMultipleScattering;
        theParticleProcessManager=theAntiProtonProcessManager;
        outFile << " ---------- particle = antiproton ----------------" << endl;
     }
     else
     {
      G4cout << " Do you want the pi+ as particle (yes/no)? " << flush;
      cin >> confirm ;
      if(confirm == "yes")
      {
      theParticle = thePionPlus;
      theParticleMultipleScattering=&thePionPlusMultipleScattering;
      theParticleProcessManager=thePionPlusProcessManager;
      outFile << " ---------- particle = pi+ ----------------" << endl;
      }
      else
      {
        G4cout << " Do you want the pi- as particle (yes/no)? " << flush;
        cin >> confirm ;
        if(confirm == "yes")
        {
        theParticle = thePionMinus;
        theParticleMultipleScattering=&thePionMinusMultipleScattering;
        theParticleProcessManager=thePionMinusProcessManager;
        outFile << " ---------- particle = pi- ----------------" << endl;
        } 
        else
        {
          G4cout << " Do you want the K+ as particle (yes/no)? " << flush;
          cin >> confirm ;
          if(confirm == "yes")
          {
          theParticle = theKaonPlus;
          theParticleMultipleScattering=&theKaonPlusMultipleScattering;
          theParticleProcessManager=theKaonPlusProcessManager;
          outFile << " ---------- particle = K+ ----------------" << endl;
          }
          else
          {
            G4cout << " Do you want the K- as particle (yes/no)? " << flush;
            cin >> confirm ;
            if(confirm == "yes")
            {
            theParticle = theKaonMinus;
            theParticleMultipleScattering=&theKaonMinusMultipleScattering;
            theParticleProcessManager=theKaonMinusProcessManager;
            outFile << " ---------- particle = K- ----------------" << endl;
            }
            else
            {
             G4cout << " There is no other particle in the test." << endl;
             return EXIT_FAILURE;
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


  G4cout << "give cuts in range" << endl ;

  G4cout << "cut for GAMMA in mm =" ;
  cin >> cutinrange ; 
  theGamma->SetCuts(cutinrange) ;
    G4cout << "gamma,cut in range(mm)=" << theGamma->GetCuts() << endl ;
    outFile << "  ---------------------------------------" << endl ;
    outFile << "  gamma,cut in range(mm)=" << theGamma->GetCuts() << endl ;

  GammaKineticEnergyCuts = theGamma->GetCutsInEnergy() ;
  for (G4int icut=0; icut<theMaterialTable->length(); icut++)
  {
    G4cout << "material index=" << icut << "  kin.energy cut(MeV)=" << 
         GammaKineticEnergyCuts[icut] << endl ;
    outFile << "  material index=" << icut << "  kin.energy cut(MeV)=" << 
         GammaKineticEnergyCuts[icut] << endl ;
  }

  G4cout << "cut for ELECTRON in mm =" ;
  cin >> cutinrange ; 
  theElectron->SetCuts(cutinrange) ;
    G4cout << "electron,cut in range(mm)=" << theElectron->GetCuts() << endl ;
    outFile << "  ---------------------------------------" << endl ;
    outFile << "  electron,cut in range(mm)=" << theElectron->GetCuts() << endl ;

  ElectronKineticEnergyCuts = theElectron->GetCutsInEnergy() ;
  for ( icut=0; icut<theMaterialTable->length(); icut++)
  {
    G4cout << "material index=" << icut << "  kin.energy cut(MeV)=" << 
         ElectronKineticEnergyCuts[icut] << endl ;
    outFile << "  material index=" << icut << "  kin.energy cut(MeV)=" << 
         ElectronKineticEnergyCuts[icut] << endl ;
  }

  G4cout << "cut for POSITRON in mm =" ;
  cin >> cutinrange ; 
  thePositron->SetCuts(cutinrange) ;
    G4cout << "positron,cut in range(mm)=" << thePositron->GetCuts() << endl ;
    outFile << "  ---------------------------------------" << endl ;
    outFile << "  positron,cut in range(mm)=" << thePositron->GetCuts() << endl ;

  PositronKineticEnergyCuts = thePositron->GetCutsInEnergy() ;
  for ( icut=0; icut<theMaterialTable->length(); icut++)
  {
    G4cout << "material index=" << icut << "  kin.energy cut(MeV)=" << 
         PositronKineticEnergyCuts[icut] << endl ;
    outFile << "  material index=" << icut << "  kin.energy cut(MeV)=" << 
         PositronKineticEnergyCuts[icut] << endl ;
  }

//****************************************************
// setcut for the selected particle (if it is not e-/e+)
 if((theParticle != theElectron) && (theParticle != thePositron))
 {
  G4cout << "cut for the selected particle in mm =" ;
  cin >> cutinrange ; 
  theParticle->SetCuts(cutinrange) ;


    G4cout << "PARTICLE: cut in range(mm)=" << theParticle->GetLengthCuts() << endl ;
    outFile << "  ---------------------------------------" << endl ;
    outFile << "PARTICLE: cut in range(mm)=" << theParticle->GetLengthCuts() << endl ;

  ParticleKineticEnergyCuts = theParticle->GetEnergyCuts() ;

  for ( icut=0; icut<theMaterialTable->length(); icut++)
  {
    G4cout << "material index=" << icut << "  kin.energy cut(MeV)=" << 
         ParticleKineticEnergyCuts[icut] << endl ;
    outFile << "  material index=" << icut << "  kin.energy cut(MeV)=" << 
         ParticleKineticEnergyCuts[icut] << endl ;
  }
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
    
    outFile << "  " << endl;
    outFile << " M S C test **********************************************" << endl ;
    outFile << "  " << endl;
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
  const G4int Nbin=97 ;
  G4double TkinMeV[Nbin]  =
              {0.001,0.0015,0.002,0.003,0.004,0.005,0.006,0.008,
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
               1.0e9} ; 
         
    G4int J=-1 ;

    G4double lambda,trueStep,geomStep,stepLimit,
             previousStepSize,currentMinimumStep ;
    G4ParticleChange* aParticleChange ;
    

    NEXTMATERIAL: ;
    J = J+1 ;
    if ( J >= theMaterialTable->length() )
      { G4cout << "that was the last material in the table --> STOP" << endl;
        return EXIT_FAILURE ; }  

    apttoMaterial = (*theMaterialTable)[ J ] ;
    MaterialName = apttoMaterial->GetName() ; 
    G4cout << "material=" << MaterialName << endl ;
    G4cout << "Do you want the MSC  Test for this material?" << endl ;
    G4cout << "type a positive number if the answer is YES" << endl ;
    G4cout << "type a negative number if the answer is NO " << endl ;
    G4int icont ;
    cin >> icont ;
    if ( icont < 0 )
        goto NEXTMATERIAL ;

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

    G4cout <<  endl;
    G4cout <<"  " << MaterialName  << "  Along Step test" << endl ;
    G4cout << " ++++++++++++++++++++++++++++++++++++++++++++" << endl ;
    G4cout << endl ;
    G4cout << "kin.en.(MeV)    lambda(mm)    trueStep(mm)" ;
    G4cout << "--->geomStep(mm)--->trueStep(mm)" << endl ;
    G4cout << endl ;
 
    outFile <<  endl;
    outFile <<"  " << MaterialName  << "  Along Step test" << endl ;
    outFile << " +++++++++++++++++++++++++++++++++++++++++" << endl ;
    outFile << endl ;
    outFile << "kin.en.(MeV)    lambda(mm)    trueStep(mm)" ;
    outFile << "--->geomStep(mm)--->trueStep(mm)" << endl ;
    outFile << endl ;
 

    for ( G4int i=0 ; i<Nbin ; i++)
    {
      trueStep = cutinrange ; 
      previousStepSize = cutinrange ;
      currentMinimumStep = trueStep ;
      (*tracke).SetKineticEnergy(TkinMeV[i]) ;
      stepLimit = (*palongget)(1)->AlongStepGetPhysicalInteractionLength( 
                                                         trackele,            
                                                         previousStepSize,
                                                         currentMinimumStep,
                                                         currentSafety,
                                                         &selection) ;

      lambda = theParticleMultipleScattering->GetTransportMeanFreePath();
      geomStep = stepLimit ;
      (*tracke).SetStepLength(geomStep) ;
      aParticleChange = (*palongdo)(0)->AlongStepDoIt(
                                               trackele,                    
                                               Step) ;

       G4cout <<" " <<  TkinMeV[i] << "  " << lambda/mm << "  " ;
       G4cout << trueStep/mm << "  " << geomStep/mm << "  " ;
       G4cout << (*aParticleChange).GetTrueStepLength()/mm << endl ;

       outFile <<" " <<  TkinMeV[i] << "  " << lambda/mm << "  " ;
       outFile << trueStep/mm << "  " << geomStep/mm << "  " ;
       outFile << (*aParticleChange).GetTrueStepLength()/mm << endl ;
    }

    G4cout <<  endl;
    outFile << endl;

   SCATTERING: ;

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

    G4cout << endl ;
    G4cout << "test of PostStepDoIt (scattering+lateral displacement) comes" ;
    G4cout << endl ;
    G4cout << "type a positive number if you want it" << endl;
    cin >> icont ;
    if ( icont<=0) goto TIMING ;

    outFile << endl ;
    outFile << "  " << MaterialName << "  PostStepDoIt (scattering) test " << endl ;
    outFile << "  +++++++++++++++++++++++++++++++++++++++++++++++++" << endl ;
    outFile << endl ;

    G4double wg[100] ;

    G4double over,overx ;
    G4double fwg, costheta,theta ;
    G4int ibin,iw ;
    G4double TMeV ;
    G4double dthetadeg,dtheta ;
    G4int flagdeg ;

   G4cout << " Do you want to give the kinetic energy (yes/no) ?" << flush;
   cin >> confirm ;
   if(confirm == "yes")
   {
    G4cout << " give the kinetic energy in MeV: " ;
    cin >> TMeV ;
   }
   else
   {
    G4cout << " give the particle momentum in MeV: " ;
    cin >> TMeV ;
    TMeV = sqrt(TMeV*TMeV+mass*mass)-mass ;
   } 

    G4cout << " give the (geom.) Step in mm: " ;
    G4double gstepmm ;
    cin >> gstepmm ;
  
    G4cout << " give number of events you want: " ;
    G4int events ;
    cin >> events ;

   G4cout << " Do you want to give the angles in degree (yes/no) ?" << flush;
   cin >> confirm ;
   if(confirm == "yes")
   {
    G4cout << " give width of the theta bin in degree:" ;
    cin >> dthetadeg ;  
    dtheta = twopi*dthetadeg/360. ;
    flagdeg = 1 ;


    G4cout << "give theta1/e from experiment (in deg):" ;
    cin >> theta1edata ;
    theta1edata*=twopi/360. ;
    G4cout << " give error of this value (in deg):" ;
    cin >> errth ;
    errth *= twopi/360. ;
   }
   else
   {
    G4cout << " give width of the theta bin in radian:" ;
    cin >> dtheta ;
    flagdeg = 0 ;
    dthetadeg = dtheta*360./twopi ;
    G4cout << " give theta1/e from experiment (in rad):" ;
    cin >> theta1edata ;
    G4cout << " give error of this value (in rad):" ;
    cin >> errth ;
   } 

//*************************************************************************************
// example: how to use the GetRange() fuctions............

   G4double Charge,rangestart,stepoverrange ;
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

   stepoverrange = gstepmm/rangestart ;
   G4cout << endl ;
   G4cout << " range at start= " << rangestart << 
           "  geom.Step/range=" <<  stepoverrange << endl ;
   G4cout << endl ;   
   outFile << endl ;
   outFile << " range at start= " << rangestart << 
           "  geom.Step/range=" <<  stepoverrange << endl ;

   outFile << endl ;   
//*************************************************************************************

    fwg=twopi/(4.*180.*180.) ;
    fwg /= events ;

    over = 0. ;
    overx = 0. ;
    thetamean = 0. ;
    thetamean2 = 0. ;

    for ( iw=0; iw<100; iw++)
    {
      theta=iw*dtheta ;
      wg[iw]=fwg/(cos(theta)-cos(theta+dtheta)) ;
    }

   //  compute TransportMeanFreePath first
       trueStep=cutinrange ;
       previousStepSize=cutinrange ;
       currentMinimumStep=trueStep ;
       (*tracke).SetKineticEnergy(TMeV) ;
       stepLimit=(*palongget)(1)->AlongStepGetPhysicalInteractionLength(
                                                       trackele,
                                                       previousStepSize,
                                                       currentMinimumStep,
                                                       currentSafety,
                                                       &selection ) ;
   // compute true path length from geom
       (*tracke).SetStepLength(gstepmm) ;
       aParticleChange = (*palongdo)(0)->AlongStepDoIt(
                                             trackele,
                                             Step) ;
       trueStep=(*aParticleChange).GetTrueStepLength();
      (*Step).SetTrack(tracke) ;
      (*Step).SetStepLength(trueStep);


   //  now you can test PostStepDoIt
       aParticleChange = (*ppostdo)(0)->PostStepDoIt(
                                             trackele,
                                             Step) ;
       const G4ThreeVector* finalPos ;
       finalPos = (*aParticleChange).GetPositionChange() ;
       xc=(*finalPos).x() ;
       yc=(*finalPos).y() ;
       zc=(*finalPos).z() ;
       later = sqrt(xc*xc+yc*yc+zc*zc) ; 


   G4cout << endl ;
   G4cout << "kin.energy=" << TMeV << " MeV   " << "geom.Step=" << gstepmm << " mm" ;
   G4cout << endl ;
   G4cout << "  lateral displacement=" << later << " mm" << endl ;
   outFile << "  mean lateral distribution -----------------" << endl ;
   outFile << endl ;
   outFile << "kin.energy=" << TMeV << " MeV   " << "geom.Step=" << gstepmm << " mm" ;
   outFile << endl ;
   outFile << "  lateral displacement=" << later << " mm" << endl ;

   
   lambda = theParticleMultipleScattering->GetTransportMeanFreePath();
   G4cout << "true Step length(mm)=" << trueStep ;
   G4cout << "   lambda(mm)=" << lambda << endl;
   outFile << "true Step length(mm)=" << trueStep ;
   outFile << "   lambda(mm)=" << lambda << endl;

  //  loop on events 

   theTimer.Start() ;

     for ( ev=0 ; ev<events; ev++)
     {
       (*tracke).SetStepLength(gstepmm) ;
       aParticleChange = (*palongdo)(0)->AlongStepDoIt(
                                             trackele,
                                             Step) ;
       trueStep=(*aParticleChange).GetTrueStepLength();
      (*Step).SetTrack(tracke) ;
      (*Step).SetStepLength(trueStep);

       aParticleChange = (*ppostdo)(0)->PostStepDoIt(
                                             trackele,
                                             Step) ;
       const G4ThreeVector* finalDir ;
       finalDir = (*aParticleChange).GetMomentumChange() ;
       costheta = (*finalDir).z() ;
       theta = acos(costheta) ;

       thetamean += theta ;
       thetamean2 += theta*theta ;

       ibin = theta/dtheta ;
       if(ibin>99)
         over += 1. ;
       else
       {
         distr[ibin] += wg[ibin] ;
         nev[ibin] += 1 ;
       }
       if((*finalDir).x()>=0.)
       {
       costheta=costheta/sqrt(costheta*costheta+((*finalDir).x())*((*finalDir).x()));
       theta = acos(costheta);
       ibin=theta/dtheta ;
       if(ibin>99)
         overx +1. ;
       else
       {
         distrx[ibin] += 1. ;
         nevx[ibin] += 1 ;
       }
       }  
     }

   theTimer.Stop() ;
 
    thetamean /= events ;
    thetamean2 /= events ;
    dthetamean = sqrt((thetamean2-thetamean*thetamean)/events) ;

     G4cout << endl ;
     G4cout << "  theta distribution  kin.energy=" << TMeV << " MeV" ;
     G4cout << "  geom.Step=" << gstepmm << " mm" << endl ;
     G4cout << " nb of events=" << events << " nb of overflows=" << over << endl;
     G4cout << " execution time/event=" << theTimer.GetUserElapsed()/events
          << " sec" << endl;
     G4cout << endl ;

     outFile << endl ;
     outFile << "  theta distribution  kin.energy=" << TMeV << " MeV" ;
     outFile << "  geom.Step=" << gstepmm << " mm" << endl ;
     outFile << " nb of events=" << events << " nb of overflows=" << over << endl;
     outFile << " execution time/event=" << theTimer.GetUserElapsed()/events 
             << " sec" << endl;
     outFile << endl ;

     theta1e=distr[0]/exp(1.) ;
     i1e=-1 ;

    G4cout << endl ;
    G4cout << "    thetalow      distr(space)                    distr(projected) "
         << endl;
    outFile << endl ;
    outFile << "    thetalow      distr(space)                    distr(projected) "
            << endl;

     for ( ib=0 ; ib<100 ; ib++)
     {
       if(flagdeg == 1)
       theta=ib*dthetadeg ;
       if(flagdeg == 0)
       theta=ib*dtheta ;

       if(distr[ib]>theta1e)
          i1e=ib ;

       errdistr = DBL_MAX ;
       if( nev[ib]>0)
          errdistr = distr[ib]/sqrt(nev[ib]) ;

       G4cout << "  " << theta << "   " << distr[ib] << 
               " +- " << errdistr ;
       outFile << "  " << theta << "   " << distr[ib] <<
               " +- " << errdistr ;
       errdistr = DBL_MAX ;
       if(nevx[ib]>0)
          errdistr = distrx[ib]/sqrt(nevx[ib]) ;
       G4cout << "    " << distrx[ib] << " +- " << errdistr << endl;
       outFile << "    " << distrx[ib] << " +- " << errdistr << endl;

     }
     if(i1e>-1)
     {
      G4cout << endl ;
      outFile << endl ;
      G4cout << "mean scattering angle in radian = " <<
              thetamean << " +- " << dthetamean << endl ;
      outFile << "mean scattering angle in radian = " <<
              thetamean << " +- " << dthetamean << endl ;
      thetamean = 360.*thetamean/twopi ;
      dthetamean = 360.*dthetamean/twopi ;
      G4cout << "mean scattering angle in degree = " <<
              thetamean << " +- " << dthetamean << endl ;
      outFile << "mean scattering angle in degree = " <<
              thetamean << " +- " << dthetamean << endl ;
      G4cout << endl ;
      outFile << endl ;

      G4cout << " theta(1/e) from the distribution:  " << endl;
      outFile << " theta(1/e) from the distribution:  " << endl;
       th1=i1e*dthetadeg ;
       th2=th1+dthetadeg ;
       G4cout << endl ;
       G4cout << "in deg.  " << th1 << " < theta(1/e) < " << th2 << endl ;
       outFile << endl ;
       outFile << "in deg.  " << th1 << " < theta(1/e) < " << th2 << endl ;
       dtheta=dthetadeg*3.1415927/180. ;
       th1=i1e*dtheta ;
       th2=th1+dtheta ;
       G4cout << "in rad.  " << th1 << " < theta(1/e) < " << th2 << endl ;
       outFile << "in rad.  " << th1 << " < theta(1/e) < " << th2 << endl; 
       G4cout << endl;
       outFile << endl;
       th1 = sqrt(2.*(exp(1./3.)-1.)*trueStep/lambda) ;
       th2 = th1*360./twopi ;
       G4cout << " theta(1/e) from the model function (theta<<1.appr.) " << endl;
       G4cout << " in degree: " << th2 << "     in radian:" << th1 << endl;
       outFile << "theta(1/e) from the model function (theta<<1.appr.)" << endl;
       outFile << " in degree: " << th2 << "     in radian:" << th1 << endl;
       th1 = acos(1.-(exp(1./3.)-1.)*trueStep/lambda) ;
       th2 = th1*360./twopi ;
       G4cout << " theta(1/e) from the model function  " << endl;
       G4cout << " in degree: " << th2 << "     in radian:" << th1 << endl;
       outFile << "theta(1/e) from the model function " << endl;
       outFile << " in degree: " << th2 << "     in radian:" << th1 << endl;
     }
     G4cout << endl ;
     outFile << endl ;

     lambdadata=(exp(1./3.)-1.)*trueStep/(1.-cos(theta1edata)) ;
     lambdadatamax=(exp(1./3.)-1.)*trueStep/(1.-cos(theta1edata-errth)) ;
     lambdadatamin=(exp(1./3.)-1.)*trueStep/(1.-cos(theta1edata+errth)) ;
     G4cout << endl;
     G4cout << " lambda from the experimental theta(1/e) (in mm) " << endl;
     G4cout << "     min.            mean              max." << endl;
     G4cout << lambdadatamin << "   " << lambdadata << "   " << lambdadatamax << endl;
     G4cout << endl ;
     outFile << endl;
     outFile << " lambda from the experimental theta(1/e) (in mm) " << endl;
     outFile << "     min.            mean              max." << endl;
     outFile << lambdadatamin << "   " << lambdadata << "   " << lambdadatamax << endl;
     outFile << endl ;
            
   goto SCATTERING ;

   TIMING: ;
 // timing test ............................................

    G4cout << " timing test follows .............." << endl ;
    G4cout << endl ;
    G4cout << "type a positive number if you want it" << endl;
    cin >> icont ;
    if ( icont<=0) goto NEXTMATERIAL ;

    outFile << endl ;
    outFile << "  " << MaterialName << "  TIMING test  " << endl ;
    outFile << "  +++++++++++++++++++++++++++++++++++++++++++++++++" << endl ;
    outFile << endl ;

    G4cout << " give the kinetic energy in MeV: " ;
    cin >> TMeV ;

    G4cout << " give the (geom.) Step in mm: " ;
    cin >> gstepmm ;
  
    G4cout << " give number of events you want: " ;
    cin >> events ;

   theTimer.Start() ;

   for ( ev=0; ev<events; ev++)
   {

   //  compute TransportMeanFreePath first
       trueStep=cutinrange ;
       previousStepSize=cutinrange ;
       currentMinimumStep=trueStep ;
       (*tracke).SetKineticEnergy(TMeV) ;
       stepLimit=(*palongget)(1)->AlongStepGetPhysicalInteractionLength(
                                                       trackele,
                                                       previousStepSize,
                                                       currentMinimumStep,
                                                       currentSafety,
                                                       &selection) ;
 
   // compute true path length from geom
       (*tracke).SetStepLength(gstepmm) ;
       aParticleChange = (*palongdo)(0)->AlongStepDoIt(
                                             trackele,
                                             Step) ;
 
   // call PostStepGetPhysicalInteractionLength
       stepLimit=(*ppostget)(1)->PostStepGetPhysicalInteractionLength(
                                                        trackele,
                                                        previousStepSize,
                                                        condition) ;
 
       trueStep=(*aParticleChange).GetTrueStepLength();
      (*Step).SetTrack(tracke) ;
      (*Step).SetStepLength(trueStep);


   //  call PostStepDoIt
       aParticleChange = (*ppostdo)(0)->PostStepDoIt(
                                             trackele,
                                             Step) ;
 
            
   }
 
   theTimer.Stop() ;

     G4cout << endl ;
     G4cout << "  timing  test  kin.energy=" << TMeV << " MeV" ;
     G4cout << "  geom.Step=" << gstepmm << " mm" << endl ;
     G4cout << " nb of events=" << events <<  endl;
     G4cout << " execution time=/event" << theTimer.GetUserElapsed()/events
          << " sec" << endl;
     G4cout << endl ;

     outFile << endl ;
     outFile << "  timing test  kin.energy=" << TMeV << " MeV" ;
     outFile << "  geom.Step=" << gstepmm << " mm" << endl ;
     outFile << " nb of events=" << events <<  endl;
     outFile << endl ;
     outFile << " AlongStepGetPhysicalInteractionLength " << endl ;
     outFile << " AlongStepDoIt " << endl ;
     outFile << " PostStepGetPhysicalInteractionLength " << endl ;
     outFile << " PostStepDoIt have been called for every events! " << endl ;
     outFile << endl ;
     outFile << " execution time/event=" << theTimer.GetUserElapsed()/events
             << " sec" << endl;
     outFile << endl ;

  
   goto TIMING ;

    if( J < theMaterialTable->length()-1 )
       goto NEXTMATERIAL ;
 
  return EXIT_SUCCESS;
}
