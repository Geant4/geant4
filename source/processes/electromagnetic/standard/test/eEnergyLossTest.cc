// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: eEnergyLossTest.cc,v 1.1 1999-01-08 16:32:26 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//-----------------------------------------------------------------
#include "G4ios.hh"
#include <fstream.h>
#include <iomanip.h>
#include "g4templates.hh"
#include "globals.hh"
#include "G4Timer.hh"     
#include "G4eEnergyLoss.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4GRSVolume.hh"
#include "G4Box.hh"
#include "G4ProcessManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4GPILSelection.hh"

//   It tests the G4eEnergyLoss,G4eIonisation,G4eBremsstrahlung processes -----------
//    created by L.Urban on 21/03/97 --------------------------
//
//    Modifications:
//    23-09-97: geometry adapted for the touchable.
//    22-04-97: adapted for G4VparticleChange, ProcessVectorTypeIndex.  mma
//

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
  //---write results to the file  ion.out-----
   ofstream outFile("eloss", ios::out ) ;
   outFile.setf( ios::scientific, ios::floatfield );

  //--------- Material definition ---------
  G4double a, z, ez, density ,temperature,pressure;
  G4State state ;
  G4String name, symbol;
  G4int nel;

  a = 9.012*g/mole;
  density = 1.848*g/cm3;
  G4Material* Be = new G4Material(name="Beryllium", z=4. , a, density);

  a = 26.98*g/mole;
  density = 2.7*g/cm3;
  G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);

  a = 28.09*g/mole;
  density = 2.33*g/cm3;
  G4Material* Si = new G4Material(name="Silicon", z=14., a, density);

  G4Element*   elH = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);

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

  G4Material* H2O = new G4Material ("Water" , 1.*g/cm3, 2);
  H2O->AddElement(elH,2);
  H2O->AddElement(elO,1);

  a = 55.85*g/mole;
  density = 7.87*g/cm3;
  G4Material* Fe = new G4Material(name="Iron", z=26., a, density);

  a = 196.97*g/mole;
  density = 19.32*g/cm3;
  G4Material* Au = new G4Material(name="Gold", z=79., a, density);

  a = 207.19*g/mole;
  density = 11.35*g/cm3;
  G4Material* Pb = new G4Material(name="Lead", z=82., a, density);

  //VacuumOnEarth  (air with a vely low density , like in the beam pipe of an accelerator)
  density = 1.0e-10*g/cm3;
  G4Material* VacuumOnEarth = new G4Material(name="VacuumOnEarth", density, nel=2);
  VacuumOnEarth->AddElement(elN, .7);
  VacuumOnEarth->AddElement(elO, .3);

  //VacuumInSpace  (H , density is extremely small)
  a = 1.01*g/mole;
  density = 1.0e-50*g/cm3;
  G4Material* VacuumInSpace = new G4Material(name="VacuumInSpace",
                            z=1., a, density) ;


  const G4MaterialTable* theMaterialTable ;
  G4Material* apttoMaterial ;
  G4String MaterialName ; 
  G4Timer theTimer ;
//--------- Particle definition ---------

  G4Gamma* theGamma = G4Gamma::GammaDefinition();
  G4Electron* theElectron = G4Electron::ElectronDefinition();
  G4Positron* thePositron = G4Positron::PositronDefinition();

  G4double* ElectronKineticEnergyCuts ;
  G4double* PositronKineticEnergyCuts ;
  G4double* GammaKineticEnergyCuts ;
  G4double* ParticleKineticEnergyCuts ;

  theMaterialTable = G4Material::GetMaterialTable() ;

  G4double cutinrange,CutInRangeele,CutInRangepos ;

  G4eIonisation theElectronIonisation , thePositronIonisation ;
  G4ProcessManager* theElectronProcessManager = theElectron->GetProcessManager();
  theElectronProcessManager->AddProcess(&theElectronIonisation,-1,0,0) ;
  G4ProcessManager* thePositronProcessManager = thePositron->GetProcessManager();
  thePositronProcessManager->AddProcess(&thePositronIonisation,-1,0,0) ;

  G4eBremsstrahlung theElectronBremsstrahlung , thePositronBremsstrahlung ;
  theElectronProcessManager->AddProcess(&theElectronBremsstrahlung,-1,-1,1) ;
  thePositronProcessManager->AddProcess(&thePositronBremsstrahlung,-1,-1,1) ;

  G4ForceCondition cond ;
  G4ForceCondition* condition = &cond ;

  G4GPILSelection selection;
  G4double currentSafety ;
  G4double& refsafety=currentSafety;


  G4cout << "cut for GAMMA in mm =" ;
  cin >> cutinrange ; cutinrange *= mm;
  theGamma->SetCuts(cutinrange) ;
    G4cout << "gamma,cut in range(mm)=" << theGamma->GetCuts()/mm << endl ;
    outFile << "  ---------------------------------------" << endl ;
    outFile << "  gamma,cut in range(mm)=" << theGamma->GetCuts()/mm << endl ;

  GammaKineticEnergyCuts = theGamma->GetCutsInEnergy() ;
  for (G4int icut=0; icut<theMaterialTable->length(); icut++)
  {
    G4cout << "material index=" << icut << "  kin.energy cut(MeV)=" << 
         GammaKineticEnergyCuts[icut]/MeV << endl ;
    outFile << "  material index=" << icut << "  kin.energy cut(MeV)=" << 
         GammaKineticEnergyCuts[icut]/MeV << endl ;
  }

  G4cout << "cut for ELECTRON in mm =" ;
  cin >> cutinrange ; cutinrange *= mm;
  CutInRangeele = cutinrange ;
  theElectron->SetCuts(cutinrange) ;
    G4cout << "electron,cut in range(mm)=" << theElectron->GetCuts()/mm << endl ;
    outFile << "  ---------------------------------------" << endl ;
    outFile << "  electron,cut in range(mm)=" << theElectron->GetCuts()/mm << endl ;

  ElectronKineticEnergyCuts = theElectron->GetCutsInEnergy() ;
  for ( icut=0; icut<theMaterialTable->length(); icut++)
  {
    G4cout << "material index=" << icut << "  kin.energy cut(MeV)=" << 
         ElectronKineticEnergyCuts[icut]/MeV << endl ;
    outFile << "  material index=" << icut << "  kin.energy cut(MeV)=" << 
         ElectronKineticEnergyCuts[icut]/MeV << endl ;
  }

  G4cout << "cut for POSITRON in mm =" ;
  cin >> cutinrange ; cutinrange *= mm;
  CutInRangepos = cutinrange;
  thePositron->SetCuts(cutinrange) ;
    G4cout << "positron,cut in range(mm)=" << thePositron->GetCuts()/mm << endl ;
    outFile << "  ---------------------------------------" << endl ;
    outFile << "  positron,cut in range(mm)=" << thePositron->GetCuts()/mm << endl ;

  PositronKineticEnergyCuts = thePositron->GetCutsInEnergy() ;
  for ( icut=0; icut<theMaterialTable->length(); icut++)
  {
    G4cout << "material index=" << icut << "  kin.energy cut(MeV)=" << 
         PositronKineticEnergyCuts[icut]/MeV << endl ;
    outFile << "  material index=" << icut << "  kin.energy cut(MeV)=" << 
         PositronKineticEnergyCuts[icut]/MeV << endl ;
  }


    
  G4double energy, momentum, mass;
  G4ProcessVector* palongget ;
  G4ProcessVector* palongdo ;
  G4ProcessVector* ppostget ;
  G4ProcessVector* ppostdo ;

    mass   = theElectron->GetPDGMass();
    energy = 1.*GeV + mass ;
    momentum=sqrt(energy*energy-mass*mass) ;
  
    G4ParticleMomentum theMomentum(momentum,0.,0.);
  
    G4double pModule = theMomentum.mag();

    G4DynamicParticle anElectron(theElectron,energy,theMomentum);

    anElectron.SetKineticEnergy(energy-mass);

    G4DynamicParticle aPositron(thePositron,energy,theMomentum);

    aPositron.SetKineticEnergy(energy-mass);

// -------------------------------------------------------------------
  G4cout << " type a positive number if you want the test with POSITRON " << endl ;
  G4cout << " type a negative number if you want the test with ELECTRON " << endl ;
  G4double charge ;

  cin >> charge ;

  if(charge<0.)
  {
    G4cout << "  ------   E L E C T R O N  ----- " << endl ;
    outFile << "  " << endl;
    outFile << "     ELECTRON ionisation test  **************************************" << endl ;
    outFile << "  " << endl;
    ParticleKineticEnergyCuts = ElectronKineticEnergyCuts ;
    palongget = anElectron.GetDefinition()->GetProcessManager()
                                 ->GetAlongStepProcessVector(typeGPIL);
    ppostget = anElectron.GetDefinition()->GetProcessManager()
                                 ->GetPostStepProcessVector(typeGPIL);
    palongdo = anElectron.GetDefinition()->GetProcessManager()
                                 ->GetAlongStepProcessVector(typeDoIt);
    ppostdo = anElectron.GetDefinition()->GetProcessManager()
                                 ->GetPostStepProcessVector(typeDoIt);
  }
  else
  {       
    G4cout << "  ----- P O S I T R O N  ----- " << endl ;
    outFile << "  " << endl;
    outFile << "     POSITRON ionisation test  **************************************" << endl ;
    outFile << "  " << endl;
    ParticleKineticEnergyCuts = PositronKineticEnergyCuts ;
    palongget = aPositron.GetDefinition()->GetProcessManager()
                                 ->GetAlongStepProcessVector(typeGPIL);
    ppostget = aPositron.GetDefinition()->GetProcessManager()
                                 ->GetPostStepProcessVector(typeGPIL);
    palongdo = aPositron.GetDefinition()->GetProcessManager()
                                 ->GetAlongStepProcessVector(typeDoIt);
    ppostdo = aPositron.GetDefinition()->GetProcessManager()
                                 ->GetPostStepProcessVector(typeDoIt);
  }

//---------------------------------- Physics --------------------------------

  G4int itry=1, Ntry=1, Nstart, ir;
  G4double r ;

//**************************************************************************
  const G4int Nbin=129 ;
  G4double TkinMeV[Nbin]  =
              {1.e-7,1.5e-7,2.e-7,3.e-7,4.e-7,5.e-7,6.e-7,8.e-7,
               1.e-6,1.5e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,8.e-6,
               0.00001,0.000015,0.00002,0.00003,0.00004,0.00005,0.00006,0.00008,
               0.0001,0.00015,0.0002,0.0003,0.0004,0.0005,0.0006,0.0008, 
               0.001,0.0015,0.002,0.003,0.004,0.005,0.006,0.008,
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
               for (G4int k=0; k<Nbin; k++) TkinMeV[k] *= MeV; 
         
    G4int J=-1 ;

    G4double lambda,trueStep,geomStep,stepLimit,
             previousStepSize,currentMinimumStep ;
    G4VParticleChange* adummy;        
    G4ParticleChange* aParticleChange ;
    
    G4double T,dEdx,range,dEdxdelta,sdEdxdelta,sloss,mloss ;

    NEXTMATERIAL: ;
    J = J+1 ;
    if ( J >= theMaterialTable->length() )
      { G4cout << "that was the last material in the table --> STOP" << endl;
        return EXIT_FAILURE ; }  

    apttoMaterial = (*theMaterialTable)[ J ] ;
    MaterialName = apttoMaterial->GetName() ; 
    G4cout << "material=" << MaterialName << endl ;
    G4cout << "Do you want the Energyloss test 1. for this material?" << endl ;
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
    G4double aTime = 0. ;

    G4Track* tracke = new G4Track(&anElectron,aTime,aPosition) ;
    G4Track& trackele = (*tracke) ;
    //(*tracke).SetVolume(myVolume) ;
    G4GRSVolume* touche = new G4GRSVolume(myVolume, NULL, aPosition);   
    (*tracke).SetTouchable(touche);    
    (*tracke).SetMomentumDirection(aDirection) ;

    G4Track* trackp = new G4Track(&aPositron,aTime,aPosition) ;
    G4Track& trackpos = (*trackp) ;
    //(*trackp).SetVolume(myVolume) ;
    G4GRSVolume* touchp = new G4GRSVolume(myVolume, NULL, aPosition);   
    (*trackp).SetTouchable(touchp);        
    (*trackp).SetMomentumDirection(aDirection) ;

    G4Step* Step = new G4Step() ;
    G4Step& Step = (*Step) ;
    tracke->SetStep(Step);
    trackp->SetStep(Step);

    G4StepPoint* aPoint = new G4StepPoint();
    (*aPoint).SetPosition(aPosition) ;
    G4double safety = 10000.*cm ;
    (*aPoint).SetSafety(safety) ;

    (*Step).SetPostStepPoint(aPoint) ;
   
//**************************************************************************

    G4cout <<  endl;
    G4cout <<"  " << MaterialName  << "  Energyloss test 1." << endl ;
    G4cout << " ++++++++++++++++++++++++++++++++++++++++++++" << endl ;
    G4cout << endl ;
    G4cout << "kin.en.(MeV)    dE/dx(MeV/mm)    range(mm)    Step(mm)" << endl ;
    G4cout << endl ;
 
    outFile <<  endl;
    outFile <<"  " << MaterialName  << "  Energyloss test 1." << endl ;
    outFile << " +++++++++++++++++++++++++++++++++++++++++" << endl ;
    outFile << endl ;
    outFile << "kin.en.(MeV)    dE/dx(MeV/mm)    range(mm)    Step(mm)" << endl ;
    outFile << endl ;
 

    for ( G4int i=0 ; i<Nbin ; i++)
    {
     if(charge<0.)
     {
      trueStep = CutInRangeele ; 
      previousStepSize = CutInRangeele ;
      currentMinimumStep = trueStep ;
      (*tracke).SetKineticEnergy(TkinMeV[i]) ;
      stepLimit = (*palongget)(0)->AlongStepGetPhysicalInteractionLength( 
                                                         trackele,            
                                                         previousStepSize,
                                                         currentMinimumStep,
                                                         refsafety,
                                                         &selection) ;
 
      dEdx = theElectronIonisation.GetdEdx() ;
      range = theElectronIonisation.GetRangeNow() ;
      T = TkinMeV[i] ;
     }
     else
     {
      trueStep = CutInRangepos ; 
      previousStepSize = CutInRangepos ;
      currentMinimumStep = trueStep ;
      (*trackp).SetKineticEnergy(TkinMeV[i]) ;
      stepLimit = (*palongget)(0)->AlongStepGetPhysicalInteractionLength( 
                                                         trackpos,            
                                                         previousStepSize,
                                                         currentMinimumStep,
                                                         refsafety,
                                                         &selection) ;

      dEdx = thePositronIonisation.GetdEdx() ;
      range = thePositronIonisation.GetRangeNow() ;
      T = TkinMeV[i] ;
     }

       G4cout <<" " <<  T/MeV << "  " << dEdx/(MeV/mm) << "  " ;
       G4cout << range/mm << "  " << stepLimit/mm << endl ; 

       outFile <<" " <<  T/MeV << "  " << dEdx/(MeV/mm) << "  " ;
       outFile << range/mm << "  " << stepLimit/mm << endl ; 

    }

    G4cout <<  endl;
    outFile << endl;

    ENERGYLOSS2: ;

    G4cout << "material=" << MaterialName << endl ;
    G4cout << "Do you want the Energyloss test 2. for this material?" << endl ;
    G4cout << "type a positive number if the answer is YES" << endl ;
    G4cout << "type a negative number if the answer is NO " << endl ;
    cin >> icont ;
    if ( icont < 0 )
        goto ENERGYLOSS3 ;

    G4double TMeV,stepmm,stepmx,meanloss,lossnow ;
  

    G4cout << "give an energy value in MeV " ;
    cin >> TMeV ; TMeV *= MeV;

    if(charge<0.)
    {
      trueStep = CutInRangeele ; 
      previousStepSize = CutInRangeele ;
      currentMinimumStep = trueStep ;
      (*tracke).SetKineticEnergy(TMeV) ;
       stepmx = (*palongget)(0)->AlongStepGetPhysicalInteractionLength(
                                              trackele,
                                              previousStepSize,
                                              currentMinimumStep,
                                              refsafety,
                                              &selection);
    }
    else
    {
      trueStep = CutInRangepos ; 
      previousStepSize = CutInRangepos ;
      currentMinimumStep = trueStep ;   
      (*trackp).SetKineticEnergy(TMeV) ;
       stepmx = (*palongget)(0)->AlongStepGetPhysicalInteractionLength(
                                              trackpos,
                                              previousStepSize,
                                              currentMinimumStep,
                                              refsafety,
                                              &selection );

    }

    G4cout << " give a steplength in mm , the max. meaningful Step is " << stepmx/mm << " mm" <<endl;
    G4cout << "Step:" ;
    cin >> stepmm ; stepmm *= mm;
 

    if(charge<0.)
    {
     (*Step).SetTrack(tracke) ;
     (*Step).SetStepLength(stepmm);
      adummy          = (*palongdo)(0)->AlongStepDoIt(trackele,Step);
      aParticleChange = (G4ParticleChange*)adummy;
      meanloss = theElectronIonisation.GetMeanLoss() ;
      lossnow = TMeV-(*aParticleChange).GetEnergyChange();
    }
    else
    {
     (*Step).SetTrack(trackp) ;
     (*Step).SetStepLength(stepmm);
      adummy          = (*palongdo)(0)->AlongStepDoIt(trackpos,Step);
      aParticleChange = (G4ParticleChange*)adummy;
      meanloss = thePositronIonisation.GetMeanLoss() ;
      lossnow = TMeV-(*aParticleChange).GetEnergyChange();
    }

    G4cout <<  endl;
    G4cout <<"  " << MaterialName  << "  Energyloss test 2." << endl ;
    G4cout << " ++++++++++++++++++++++++++++++++++++++++++++" << endl ;
    G4cout << endl ;
    G4cout << "kin.en.(MeV)    Step(mm)   meanloss(MeV)  act.loss(MeV)" << endl ;
    G4cout << TMeV/MeV << "   " << stepmm/mm << "  " << meanloss/MeV << "  " << lossnow/MeV << endl ;
    G4cout << " status change:" << (*aParticleChange).GetStatusChange() << endl ;
    G4cout << endl ;
 
    outFile <<  endl;
    outFile <<"  " << MaterialName  << "  Energyloss test 2." << endl ;
    outFile << " +++++++++++++++++++++++++++++++++++++++++" << endl ;
    outFile << endl ;
    outFile << "kin.en.(MeV)    Step(mm)   meanloss(MeV)  act.loss(MeV)" << endl ;
    outFile << TMeV/MeV << "   " << stepmm/mm << "  " << meanloss/MeV << "  " << lossnow/MeV << endl ;
    outFile << " status change:" << (*aParticleChange).GetStatusChange() << endl ;
    outFile << endl ;
 
    goto ENERGYLOSS2 ;

    ENERGYLOSS3: ;

    G4cout << "material=" << MaterialName << endl ;
    G4cout << "Do you want the Energyloss test 3. for this material?" << endl ;
    G4cout << "type a positive number if the answer is YES" << endl ;
    G4cout << "type a negative number if the answer is NO " << endl ;
    cin >> icont ;
    if ( icont < 0 )
        goto DELTARAY1 ;

    G4cout << "give an energy value in MeV " ;
    cin >> TMeV ; TMeV *= MeV;

    if(charge<0.)
    {
      trueStep = CutInRangeele ; 
      previousStepSize = CutInRangeele ;
      currentMinimumStep = trueStep ;
      (*tracke).SetKineticEnergy(TMeV) ;
       stepmx = (*palongget)(0)->AlongStepGetPhysicalInteractionLength(
                                              trackele,
                                              previousStepSize,
                                              currentMinimumStep,
                                              refsafety,
                                              &selection);
    }
    else
    {
      trueStep = CutInRangepos ; 
      previousStepSize = CutInRangepos ;
      currentMinimumStep = trueStep ;   
      (*trackp).SetKineticEnergy(TMeV) ;
       stepmx = (*palongget)(0)->AlongStepGetPhysicalInteractionLength(
                                              trackpos,
                                              previousStepSize,
                                              currentMinimumStep,
                                              refsafety,
                                              &selection);

    }

    G4cout << " give a steplength in mm , the max. meaningful Step is " << stepmx/mm << " mm" <<endl;
    G4cout << "Step:" ;
    cin >> stepmm ; stepmm *= mm;

    G4cout << " give number of events you want " ;
    G4int nbev,ibev ;
    cin >> nbev ;

    meanloss=0.;
    theTimer.Start();

    for ( ibev=0; ibev<nbev; ibev++)
    { 
    if(charge<0.)
   {
     (*Step).SetTrack(tracke) ;
     (*Step).SetStepLength(stepmm);
      adummy          = (*palongdo)(0)->AlongStepDoIt(trackele,Step);
      aParticleChange = (G4ParticleChange*)adummy;
      lossnow = TMeV-(*aParticleChange).GetEnergyChange();
    }
    else
    {
     (*Step).SetTrack(trackp) ;
     (*Step).SetStepLength(stepmm);
      adummy          = (*palongdo)(0)->AlongStepDoIt(trackpos,Step);
      aParticleChange = (G4ParticleChange*)adummy;
      lossnow = TMeV-(*aParticleChange).GetEnergyChange();
    }

    meanloss += lossnow ;
   }

    theTimer.Stop();
    meanloss /= nbev ;
    G4cout <<  endl;
    G4cout <<"  " << MaterialName  << "  Energyloss test 3." << endl ;
    G4cout << " ++++++++++++++++++++++++++++++++++++++++++++" << endl ;
    G4cout << endl ;
    G4cout << "kin.en.(MeV)    Step(mm)   meanloss(MeV) time/event(sec) " << endl ;
    G4cout << TMeV/MeV << "   " << stepmm/mm << "  " << meanloss/MeV << "  " << 
    theTimer.GetUserElapsed()/nbev << endl ;
    G4cout << endl ;
 
    outFile <<  endl;
    outFile <<"  " << MaterialName  << "  Energyloss test 3." << endl ;
    outFile << " +++++++++++++++++++++++++++++++++++++++++" << endl ;
    outFile << endl ;
    outFile << "kin.en.(MeV)    Step(mm)   meanloss(MeV) time/event(sec) " << endl ;
    outFile << TMeV/MeV << "   " << stepmm/mm << "  " << meanloss/MeV << "  " << 
    theTimer.GetUserElapsed()/nbev << endl ;
    outFile << endl ;
 
    goto ENERGYLOSS3 ;

    DELTARAY1: ;
    G4cout << "material=" << MaterialName << endl ;
    G4cout << "Do you want the delta ray test 1. for this material?" << endl ;
    G4cout << "type a positive number if the answer is YES" << endl ;
    G4cout << "type a negative number if the answer is NO " << endl ;
    cin >> icont ;
    if ( icont < 0 )
        goto DELTARAY2 ;


    G4cout <<  endl;
    G4cout <<"  " << MaterialName  << "  delta ray test 1." << endl ;
    G4cout << " ++++++++++++++++++++++++++++++++++++++++++++" << endl ;
    G4cout << endl ;
    G4cout << "kin.en.(MeV)      mean free path(mm)" << endl ;
    G4cout << endl ;
 
    outFile <<  endl;
    outFile <<"  " << MaterialName  << "  delta ray test 1." << endl ;
    outFile << " +++++++++++++++++++++++++++++++++++++++++" << endl ;
    outFile << endl ;
    outFile << "kin.en.(MeV)      mean free path(mm)" << endl ;
    outFile << endl ;

    for ( i=0 ; i<Nbin ; i++)
    {
    // G4eIonisation::StartTracking() ;

     if(charge<0.)
     {
      previousStepSize = CutInRangeele ;
      (*tracke).SetKineticEnergy(TkinMeV[i]) ;
      stepLimit = theElectronIonisation.GetMeanFreePath(                                           
                                                         trackele,            
                                                         previousStepSize,
                                                         condition) ;                                           
     }
     else
     {
      previousStepSize = CutInRangepos ;
      (*trackp).SetKineticEnergy(TkinMeV[i]) ;
      stepLimit = thePositronIonisation.GetMeanFreePath(                                           
                                                         trackpos,            
                                                         previousStepSize,
                                                         condition) ;                          
     }

      T = TkinMeV[i] ;

      G4cout <<" " <<  T/MeV << "      " <<  stepLimit/mm << endl ;

      outFile <<" " <<  T/MeV << "      " << stepLimit/mm << endl ;

     // G4eIonisation::EndTracking() ;
 
    }

    G4cout <<  endl;
    outFile << endl;

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    DELTARAY2: ;

    G4cout << "material=" << MaterialName << endl ;
    G4cout << "Do you want the deltaray test 2. for this material?" << endl ;
    G4cout << "type a positive number if the answer is YES" << endl ;
    G4cout << "type a negative number if the answer is NO " << endl ;
    cin >> icont ;

    G4double newenergy,dx,dy,dz,Tdelta,ddx,ddy,ddz ;
    G4int nd ;
    const G4ThreeVector* momdir ;
    G4ParticleMomentum ddir ;

    if ( icont < 0 )
        goto BREMS1 ;  

    G4cout << "give an energy value in MeV " ;
    cin >> TMeV ; TMeV *= MeV;
                           
     stepmm = 1.*mm ;

    if(charge<0.)
    {
     (*Step).SetTrack(tracke) ;
     (*Step).SetStepLength(stepmm);
      (*tracke).SetKineticEnergy(TMeV) ;
      adummy          = (*ppostdo)(0)->PostStepDoIt(trackele,Step);
      aParticleChange = (G4ParticleChange*)adummy;
    }
    else
    {
     (*Step).SetTrack(trackp) ;
     (*Step).SetStepLength(stepmm);
      (*trackp).SetKineticEnergy(TMeV) ;
      adummy          = (*ppostdo)(0)->PostStepDoIt(trackpos,Step);
      aParticleChange = (G4ParticleChange*)adummy;
    }

     newenergy=(*aParticleChange).GetEnergyChange() ;
     momdir=(*aParticleChange).GetMomentumChange();
     dx = (*momdir).x();
     dy = (*momdir).y();
     dz = (*momdir).z();
     nd=aParticleChange->GetNumberOfSecondaries();
 
    if(nd>0)
    {
     Tdelta=aParticleChange->GetSecondary(0)->GetKineticEnergy();
     ddir=aParticleChange->GetSecondary(0)->
                              GetMomentumDirection();
     ddx = (ddir).x();
     ddy = (ddir).y();
     ddz = (ddir).z();
    }

    G4cout <<  endl;
    G4cout <<"  " << MaterialName  << "  delta ray test 2." << endl ;
    G4cout << " ++++++++++++++++++++++++++++++++++++++++++++" << endl ;
    G4cout << endl ;
    G4cout << "T=" << TMeV/MeV << "   newT=" << newenergy/MeV << "  (MeV)" << endl ;
    G4cout << " status change:" << (*aParticleChange).GetStatusChange() << endl ;
    if(nd>0)
    G4cout << "Tdelta=" << Tdelta/MeV << endl ;
    G4cout << "new direction:" << dx << "  " << dy << "  " << dz << endl;
    if(nd>0)
    G4cout << "delta direction:" << ddx << "  " << ddy << "  " << ddz << endl ;
    G4cout << endl ;
 
    outFile <<  endl;
    outFile <<"  " << MaterialName  << "  delta ray test 2." << endl ;
    outFile << " +++++++++++++++++++++++++++++++++++++++++" << endl ;
    outFile << endl ;
    outFile << "T=" << TMeV/MeV << "   newT=" << newenergy/MeV << "   (MeV)" << endl;
    outFile << " status change:" << (*aParticleChange).GetStatusChange() << endl ;
    if(nd>0)
    outFile << "Tdelta=" << Tdelta/MeV << endl ;
    outFile << "new direction:" << dx << "  " << dy << "  " << dz << endl;
    if(nd>0)
    outFile << "delta direction:" << ddx << "  " << ddy << "  " << ddz << endl ;
    outFile << endl ;

    (*aParticleChange).Clear();


//............................................
   nbev=100000 ;
   mloss = 0. ;
   sloss = 0. ;
   theTimer.Start();

   for (ibev=0 ; ibev<nbev ; ibev++ )
    {
     if(charge<0.)
     {
      adummy          = (*ppostdo)(0)->PostStepDoIt(trackele,Step);
      aParticleChange = (G4ParticleChange*)adummy;
     }
     else
     {
      adummy          = (*ppostdo)(0)->PostStepDoIt(trackpos,Step);
      aParticleChange = (G4ParticleChange*)adummy;
     }

      nd=aParticleChange->GetNumberOfSecondaries();
      Tdelta=0. ; 
      if(nd>0)
      {
       Tdelta=aParticleChange->GetSecondary(0)->GetKineticEnergy();
      }
      mloss += Tdelta ;
      sloss += Tdelta*Tdelta ;
     (*aParticleChange).Clear();
    }
    theTimer.Stop();
    mloss /= nbev ;
    sloss /= nbev ;
    sloss = sqrt((sloss-mloss*mloss)/nbev) ;
 
    previousStepSize = cutinrange ;
    if(charge<0.)
    stepLimit = theElectronIonisation.GetMeanFreePath(                                           
                                                      trackele,            
                                                      previousStepSize,
                                                      condition) ;
    else
    stepLimit = thePositronIonisation.GetMeanFreePath(                                           
                                                      trackele,            
                                                      previousStepSize,
                                                      condition) ;
                                           
    dEdxdelta=mloss/stepLimit ;
    sdEdxdelta = sloss/stepLimit ;

    G4cout << " mean energy loss due to delta production (in MeV)=" <<
              mloss/MeV << " +- " << sloss << endl ;
    G4cout << " dE/dx due to delta production (in MeV/mm,from " <<
             nbev << " events )=" << dEdxdelta/(MeV/mm) <<
             " +- " << sdEdxdelta <<  endl ;
    G4cout << " time/delta =" << theTimer.GetUserElapsed()/nbev << endl ;
    G4cout << endl ;

            
    outFile << " mean energy loss due to delta production (in MeV)=" <<
              mloss/MeV << " +- " << sloss << endl ;
    outFile << " dE/dx due to delta production (in MeV/mm,from " <<
             nbev << " events )=" << dEdxdelta/(MeV/mm) <<
             " +- " << sdEdxdelta <<  endl ;

    outFile << " time/event=" << theTimer.GetUserElapsed()/nbev << endl ;
    outFile << endl ;
//..................................


    goto DELTARAY2 ;

    BREMS1: ;
    G4cout << "material=" << MaterialName << endl ;
    G4cout << "Do you want the brems test 1. for this material?" << endl ;
    G4cout << "type a positive number if the answer is YES" << endl ;
    G4cout << "type a negative number if the answer is NO " << endl ;
    cin >> icont ;
    if ( icont < 0 )
        goto BREMS2 ;

    G4cout <<  endl;
    G4cout <<"  " << MaterialName  << "  brems test 1." << endl ;
    G4cout << " ++++++++++++++++++++++++++++++++++++++++++++" << endl ;
    G4cout << endl ;
    G4cout << "kin.en.(MeV)      mean free path(mm)" << endl ;
    G4cout << endl ;
 
    outFile <<  endl;
    outFile <<"  " << MaterialName  << "  brems test 1." << endl ;
    outFile << " +++++++++++++++++++++++++++++++++++++++++" << endl ;
    outFile << endl ;
    outFile << "kin.en.(MeV)      mean free path(mm)" << endl ;
    outFile << endl ;

    for ( i=0 ; i<Nbin ; i++)
    {
    // G4eBremsstrahlung::StartTracking() ;

     if(charge<0.)
     {
      previousStepSize = CutInRangeele ;
      (*tracke).SetKineticEnergy(TkinMeV[i]) ;
      stepLimit = theElectronBremsstrahlung.GetMeanFreePath( 
                                                         trackele,            
                                                         previousStepSize,
                                                         condition) ;
     }
     else
     {
      previousStepSize = CutInRangepos ;
      (*trackp).SetKineticEnergy(TkinMeV[i]) ;
      stepLimit = thePositronBremsstrahlung.GetMeanFreePath( 
                                                         trackele,            
                                                         previousStepSize,
                                                         condition) ;

     }

      T = TkinMeV[i] ;

      G4cout <<" " <<  T/MeV << "      " << stepLimit/mm << endl ;

      outFile <<" " <<  T/MeV << "      " << stepLimit/mm << endl ;

    //  G4eBremsstrahlung::EndTracking() ;

    }

    G4cout <<  endl;
    outFile << endl;

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    BREMS2: ;

    G4cout << "material=" << MaterialName << endl ;
    G4cout << "Do you want the brems test 2. for this material?" << endl ;
    G4cout << "type a positive number if the answer is YES" << endl ;
    G4cout << "type a negative number if the answer is NO " << endl ;
    cin >> icont ;
    if ( icont < 0 )
        goto NEXTMATERIAL ;  

    G4cout << "give an energy value in MeV " ;
    cin >> TMeV ; TMeV *= MeV;
                           
     stepmm = 1.*mm ;

    if(charge<0.)
    {
     (*Step).SetTrack(tracke) ;
     (*Step).SetStepLength(stepmm);
      (*tracke).SetKineticEnergy(TMeV) ;
      adummy          = (*ppostdo)(1)->PostStepDoIt(trackele,Step);
      aParticleChange = (G4ParticleChange*)adummy;
    }
    else
    {
     (*Step).SetTrack(trackp) ;
     (*Step).SetStepLength(stepmm);
      (*trackp).SetKineticEnergy(TMeV) ;
      adummy          = (*ppostdo)(1)->PostStepDoIt(trackpos,Step);
      aParticleChange = (G4ParticleChange*)adummy;
    }

     newenergy=(*aParticleChange).GetEnergyChange() ;
     momdir=(*aParticleChange).GetMomentumChange();
     dx = (*momdir).x();
     dy = (*momdir).y();
     dz = (*momdir).z();
     nd=aParticleChange->GetNumberOfSecondaries();
 
    if(nd>0)
    {
     Tdelta=aParticleChange->GetSecondary(0)->GetKineticEnergy();
     ddir=aParticleChange->GetSecondary(0)->
                              GetMomentumDirection();
     ddx = (ddir).x();
     ddy = (ddir).y();
     ddz = (ddir).z();
    }

    G4cout <<  endl;
    G4cout <<"  " << MaterialName  << "  brems test 2." << endl ;
    G4cout << " ++++++++++++++++++++++++++++++++++++++++++++" << endl ;
    G4cout << endl ;
    G4cout << "T=" << TMeV/MeV << "   newT=" << newenergy/MeV << "  (MeV)" << endl ;
    G4cout << " status change:" << (*aParticleChange).GetStatusChange() << endl ;
    if(nd>0)
    G4cout << "Tgamma=" << Tdelta/MeV << endl ;
    G4cout << "new direction:" << dx << "  " << dy << "  " << dz << endl;
    if(nd>0)
    G4cout << "gamma direction:" << ddx << "  " << ddy << "  " << ddz << endl ;
    G4cout << endl ;
 
    outFile <<  endl;
    outFile <<"  " << MaterialName  << "  brems test 2." << endl ;
    outFile << " +++++++++++++++++++++++++++++++++++++++++" << endl ;
    outFile << endl ;
    outFile << "T=" << TMeV/MeV << "   newT=" << newenergy/MeV << "   (MeV)" << endl;
    outFile << " status change:" << (*aParticleChange).GetStatusChange() << endl ;

    if(nd>0)
    outFile << "Tgamma=" << Tdelta/MeV << endl ;
    outFile << "new direction:" << dx << "  " << dy << "  " << dz << endl;
    if(nd>0)
    outFile << "gamma direction:" << ddx << "  " << ddy << "  " << ddz << endl ;
    outFile << endl ;

    (*aParticleChange).Clear();


//............................................
   nbev=100000 ;
   mloss = 0. ;
   sloss = 0. ;
   theTimer.Start();

   for (ibev=0 ; ibev<nbev ; ibev++ )
    {
     if(charge<0.)
     {
      adummy          = (*ppostdo)(1)->PostStepDoIt(trackele,Step);
      aParticleChange = (G4ParticleChange*)adummy;
     }
     else
     {
      adummy          = (*ppostdo)(1)->PostStepDoIt(trackpos,Step);
      aParticleChange = (G4ParticleChange*)adummy;
     }

      nd=aParticleChange->GetNumberOfSecondaries();
      Tdelta=0. ; 
      if(nd>0)
      {
       Tdelta=aParticleChange->GetSecondary(0)->GetKineticEnergy();
      }
      mloss += Tdelta ;
      sloss += Tdelta*Tdelta ;
     (*aParticleChange).Clear();
    }
    theTimer.Stop();
    mloss /= nbev ;
    sloss /= nbev ;
    sloss = sqrt((sloss-mloss*mloss)/nbev) ;
 
    previousStepSize = cutinrange ;
    if(charge<0.)
    stepLimit = theElectronBremsstrahlung.GetMeanFreePath(                                           
                                                      trackele,            
                                                      previousStepSize,
                                                      condition) ;
    else
    stepLimit = thePositronBremsstrahlung.GetMeanFreePath(                                           
                                                      trackele,            
                                                      previousStepSize,
                                                      condition) ;
                                           
    dEdxdelta=mloss/stepLimit ;
    sdEdxdelta=sloss/stepLimit ;

    G4cout << " mean energy loss due to bremsstrahlung   (in MeV)=" <<
              mloss/MeV << " +- " << sloss << endl ;
    G4cout << " dE/dx due to bremsstrahlung   (in MeV/mm,from " <<
             nbev << " events )=" << dEdxdelta/(MeV/mm) <<
             " +- " << sdEdxdelta <<  endl ;
    G4cout << " time/delta =" << theTimer.GetUserElapsed()/nbev << endl ;
    G4cout << endl ;

            
    outFile << " mean energy loss due to bremsstrahlung   (in MeV)=" <<
              mloss/MeV << " +- " << sloss << endl ;
    outFile << " dE/dx due to bremsstrahlung   (in MeV/mm,from " <<
             nbev << " events )=" << dEdxdelta/(MeV/mm) <<
             " +- " << sdEdxdelta <<  endl ;

    outFile << " time/event=" << theTimer.GetUserElapsed()/nbev << endl ;
    outFile << endl ;
//..................................
    goto BREMS2 ;

    if( J < theMaterialTable->length()-1 )
       goto NEXTMATERIAL ;
 
  return EXIT_SUCCESS;
}

