// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyTest.cc,v 1.1 1999-06-24 19:29:41 aforti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// KaonMinusAtRestTest.cc 
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      File name:     G4LowEnergyTest
//
//      Author:        Christian Voelcker (from M. Maire)
// 
//      Creation date: ?
//
//      Modifications: 
//
//      Alessandra Forti march 1999  added ntuples
//      Alessandra Forti march 1999  added processes 
//      Alessandra Forti march 1999  added cut in energy 
// -------------------------------------------------------------------

#include "G4ios.hh"
#include <fstream.h>
#include <iomanip.h>

#include "G4Material.hh"

#include "G4ProcessManager.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyRayleigh.hh"
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4LowEnergyIonisation.hh"
#include "G4eIonisation.hh"
#include "G4hIonisation.hh"

#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"

#include "G4Step.hh"
#include "G4GRSVolume.hh"

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"

HepTupleManager* hbookManager;

main()
{

  // Setup

  G4int niter=1e4;
  G4int imat=2;
  G4int verboseLevel=1;
  G4int processID=6;

  G4cout << "How many interactions? [10], Which material? [3], which Verbose Level? [1]" << endl;
  cin >> niter >> imat >> verboseLevel;

  G4cout<<"which process?"<<endl<<setw(60)<<"[1] = G4LowEnergyPhotoElectric, [2] = G4LowEnergyCompton"<<endl;
  G4cout<<setw(60)<<"[3] = G4LowEnergyRayleigh, [4] = G4LowEnergyGammaconversion"<<endl;
  G4cout<<setw(60)<<"[5] = G4LowEnergyBremstrahlung"<<"[6] = G4LowEnergyIonisation"<<endl;

  cin >> processID;

  G4double InitEnergy = 1e-3, InitX = 0., InitY = 0., InitZ = 1.;
  G4cout<<"Enter the initial particle energy E and its direction"<<endl; 
  cin >> InitEnergy >> InitX >> InitY >> InitZ;

  G4double gammaCut = 1e-3, electronCut = 1e-3;
  G4cout<<"Set photons 1e-3 mm and electrons cuts 1e-3 mm"<<endl; 
  cin>>gammaCut>>electronCut;

  //-------- write results onto a file --------
 
  //  ofstream outFile1( "lowenergypri.out", ios::out);
  //  ofstream outFile2( "lowenergysec.out", ios::out);
  //  ofstream outFile3( "lowenergymfp.out", ios::out);

  //  outFile1.setf( ios::scientific, ios::floatfield);
  //  outFile2.setf( ios::scientific, ios::floatfield);
  //  outFile3.setf( ios::scientific, ios::floatfield);

  //  outFile1.setf(ios::left);
  //  outFile2.setf(ios::left);
  //  outFile3.setf(ios::left);

  G4cout.setf( ios::scientific, ios::floatfield );
  // -------------------------------------------------------------------

  // ALE ---- HBOOK initialization
  if(processID == 1){

  hbookManager = new HBookFile("photoelectric.hbook", 58);
  assert (hbookManager != 0);
  }
  else if(processID == 2){

    hbookManager = new HBookFile("compton.hbook", 58);
    assert (hbookManager != 0);
  }
  else if(processID == 3){

    hbookManager = new HBookFile("rayleigh.hbook", 58);
    assert (hbookManager != 0);
  }
  else if(processID == 4){

    hbookManager = new HBookFile("gammaconv.hbook", 58);
    assert (hbookManager != 0);
  }
  else if(processID == 5){

    hbookManager = new HBookFile("bremstrahlung.hbook", 58);
    assert (hbookManager != 0);
  }
  else if(processID == 6){

    hbookManager = new HBookFile("ionisation.hbook", 58);
    assert (hbookManager != 0);
  }
  // ALE ---- Book a histogram and ntuples
  G4cout<<"Hbook file name: "<<((HBookFile*) hbookManager)->filename()<<endl;
  
  // ---- primary ntuple ------
  HepTuple* ntuple1 = hbookManager->ntuple("Primary Ntuple");
  assert (ntuple1 != 0);
  
  // ---- secondaries ntuple ------
  HepTuple* ntuple2 = hbookManager->ntuple("Secondaries Ntuple");
  assert (ntuple2 != 0);
 
  // ---- mfp ntuple ------
  HepTuple* ntuple3 = hbookManager->ntuple("MeanFreePath Ntuple");
  assert (ntuple3 != 0);
 
  // ---- secondaries histos ----
  HepHistogram* hEKin;
  hEKin = hbookManager->histogram("Kinetic Energy", 100,0.,200.);
  assert (hEKin != 0);  
  
  HepHistogram* hP;
  hP = hbookManager->histogram("Momentum", 100,0.,1000.);
  assert (hP != 0);  
  
  HepHistogram* hNSec;
  hNSec = hbookManager->histogram("Number of secondaries", 40,0.,40.);
  assert (hNSec != 0);  
  
  HepHistogram* hDebug;
  hDebug = hbookManager->histogram("Debug", 100,0.,200.);
  assert (hDebug != 0);  
  

  //--------- Materials definition ---------

  G4Material* Be = new G4Material("Beryllium",    4.,  9.01*g/mole, 1.848*g/cm3);
  G4Material* Graphite = new G4Material("Graphite",6., 12.00*g/mole, 2.265*g/cm3 );
  G4Material* Al  = new G4Material("Aluminium", 13., 26.98*g/mole, 2.7 *g/cm3);
  G4Material* LAr = new G4Material("LArgon",   18., 39.95*g/mole, 1.393*g/cm3);
  G4Material* Fe  = new G4Material("Iron",      26., 55.85*g/mole, 7.87*g/cm3);
  G4Material* Cu  = new G4Material("Copper",    29., 63.55*g/mole, 8.96*g/cm3);
  G4Material*  W  = new G4Material("Tungsten", 74., 183.85*g/mole, 19.30*g/cm3);
  G4Material* Pb  = new G4Material("Lead",      82., 207.19*g/mole, 11.35*g/cm3);
  G4Material*  U  = new G4Material("Uranium", 92., 238.03*g/mole, 18.95*g/cm3);

  G4Element*   H  = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   O  = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
  G4Element*   C  = new G4Element ("Carbon"  , "C", 6. , 12.00*g/mole);
  G4Element*  Cs  = new G4Element ("Cesium"  , "Cs", 55. , 132.905*g/mole);
  G4Element*   I  = new G4Element ("Iodide"  , "I", 53. , 126.9044*g/mole);

  G4Material*  maO = new G4Material("Oxygen", 8., 16.00*g/mole, 1.1*g/cm3);

  G4Material* water = new G4Material ("Water" , 1.*g/cm3, 2);
  water->AddElement(H,2);
  water->AddElement(O,1);

  G4Material* ethane = new G4Material ("Ethane" , 0.4241*g/cm3, 2);
  ethane->AddElement(H,6);
  ethane->AddElement(C,2);
  
  G4Material* csi = new G4Material ("CsI" , 4.53*g/cm3, 2);
  csi->AddElement(Cs,1);
  csi->AddElement(I,1);

  static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  G4cout<<"The material is: "<<(*theMaterialTable)(imat)->GetName()<<endl;

  // Geometry definitions
  G4Box* theFrame = new G4Box ("Frame",92*mm, 92*mm, 92*mm);
  
  G4LogicalVolume* LogicalFrame = new G4LogicalVolume(theFrame,
						      (*theMaterialTable)(imat),
						      "LFrame", 0, 0, 0);
  
  G4PVPlacement* PhysicalFrame = new G4PVPlacement(0,G4ThreeVector(),
						   "PFrame",LogicalFrame,0,false,0);
  
  // the center-of-mass of the cube should be located at the origin!

  //--------- Particle definition ---------

  G4ParticleDefinition* electron = G4Electron::ElectronDefinition();
  G4ParticleDefinition* positron = G4Positron::PositronDefinition();
  G4ParticleDefinition* gamma = G4Gamma::GammaDefinition();
  G4ParticleDefinition* proton = G4Proton::ProtonDefinition();

  //--------- Processes definition ---------

  G4ProcessManager* theGammaProcessManager = new G4ProcessManager(gamma);
  gamma->SetProcessManager(theGammaProcessManager);

  G4ProcessManager* theElectronProcessManager = new G4ProcessManager(electron);
  electron->SetProcessManager(theElectronProcessManager);
  
  G4ProcessManager* thePositronProcessManager = new G4ProcessManager(positron);
  positron->SetProcessManager(thePositronProcessManager);
  
  G4ProcessManager* theProtonProcessManager = new G4ProcessManager(proton);
  proton->SetProcessManager(theProtonProcessManager);

  G4LowEnergyPhotoElectric  PhotoElectricProcess;
  G4LowEnergyCompton ComptonProcess;
  G4LowEnergyRayleigh RayleighProcess;
  G4LowEnergyGammaConversion GammaConvProcess;
  
  theGammaProcessManager->AddDiscreteProcess(&PhotoElectricProcess);
  theGammaProcessManager->AddDiscreteProcess(&ComptonProcess);
  theGammaProcessManager->AddDiscreteProcess(&RayleighProcess);
  theGammaProcessManager->AddDiscreteProcess(&GammaConvProcess);

  G4LowEnergyIonisation IonisationProcess;
  theElectronProcessManager->AddProcess(&IonisationProcess);
  theElectronProcessManager->SetProcessOrdering(&IonisationProcess,idxAlongStep,1);
  theElectronProcessManager->SetProcessOrdering(&IonisationProcess,idxPostStep,1);

  G4LowEnergyBremsstrahlung BremstrahlungProcess;
  theElectronProcessManager->AddProcess(&BremstrahlungProcess);
  theElectronProcessManager->SetProcessOrdering(&BremstrahlungProcess,idxAlongStep,1);
  theElectronProcessManager->SetProcessOrdering(&BremstrahlungProcess,idxPostStep,1);

  G4eIonisation IonisationPlusProcess;
  thePositronProcessManager->AddProcess(&IonisationPlusProcess);
  thePositronProcessManager->SetProcessOrdering(&IonisationPlusProcess,idxAlongStep,1);
  thePositronProcessManager->SetProcessOrdering(&IonisationPlusProcess,idxPostStep,1);

  G4hIonisation hIonisationProcess;
  theProtonProcessManager->AddProcess(&hIonisationProcess);
  theProtonProcessManager->SetProcessOrdering(&hIonisationProcess,idxAlongStep,1);
  theProtonProcessManager->SetProcessOrdering(&hIonisationProcess,idxPostStep,1);

  G4ForceCondition* condition;

  // ------- set cut and Build CrossSection Tables -------
  //
  G4Gamma::SetEnergyRange(2.5e-4*MeV,1e5*MeV);
  G4Electron::SetEnergyRange(2.5e-4*MeV,1e5*MeV);
  G4Positron::SetEnergyRange(2.5e-4*MeV,1e5*MeV);

  gamma->SetCuts(1e-6*mm);
  electron->SetCuts(1e-6*mm);
  positron->SetCuts(1e-6*mm);

  G4cout<<"the cut in energy for gamma in: "<<(*theMaterialTable)(imat)->GetName()
      <<" is: "<<G4Gamma::GetCutsInEnergy()[imat]<<endl;
  G4cout<<"the cut in energy for e- in: "<<(*theMaterialTable)(imat)->GetName()
      <<" is: "<<G4Electron::GetCutsInEnergy()[imat]<<endl;

  // -------- create 1 Dynamic Particle  ----

  G4double gammaEnergy = InitEnergy*MeV;

  G4ParticleMomentum gammaDirection(InitX,InitY,InitZ);
  
  
  G4DynamicParticle photon(G4Gamma::Gamma(),gammaDirection,gammaEnergy);
  G4DynamicParticle elecT(G4Electron::Electron(),gammaDirection,gammaEnergy);

  //--------- track definition (for this test ONLY!)------------

  G4ThreeVector aPosition(0.,0.,0.);
  G4double aTime = 0. ;

  G4Track* ptrack;
  if(processID != 5 && processID != 6)
    ptrack = new G4Track(&photon,aTime,aPosition) ;
  else
    ptrack = new G4Track(&elecT,aTime,aPosition) ;

  G4Track& aTrack = (*ptrack) ;
  
  // do I really need this?

  G4GRSVolume* touche = new G4GRSVolume(PhysicalFrame, NULL, aPosition);   
  ptrack->SetTouchable(touche);
  
  // -------- create 1 Step (for this test only)----  

  G4Step* Step = new G4Step();  
  G4Step& aStep = (*Step);
  Step->SetTrack(ptrack);
  
  // --------- check applicability
  G4ParticleDefinition* PhotonDefinition = photon.GetDefinition();
  G4ParticleDefinition* ElectronDefinition = elecT.GetDefinition();

  if(!PhotoElectricProcess.IsApplicable(*PhotonDefinition) || 
     !ComptonProcess.IsApplicable(*PhotonDefinition) || 
     !RayleighProcess.IsApplicable(*PhotonDefinition) || 
     !GammaConvProcess.IsApplicable(*PhotonDefinition) ||
     !BremstrahlungProcess.IsApplicable(*ElectronDefinition)||
     !IonisationProcess.IsApplicable(*ElectronDefinition)) {

    G4cout 
      << PhotonDefinition->GetParticleName()
      << " is not a Photon! or " 
      << ElectronDefinition->GetParticleName()
      <<" is not an Electron"<<endl;
    G4Exception("FAIL: *** exit program ***\n");
  //     return ;
  }

  //  PhotoElectricProcess.SetVerboseLevel(verboseLevel);

  // Initialize the physics tables for ALL processes
  //  IonisationProcess.MinusNbOfProcesses();
  IonisationProcess.BuildPhysicsTable(*electron);
  BremstrahlungProcess.BuildPhysicsTable(*electron);

  IonisationPlusProcess.MinusNbOfProcesses();
  IonisationPlusProcess.BuildPhysicsTable(*positron);

  hIonisationProcess.BuildPhysicsTable(*proton);

  if(processID == 1){

    PhotoElectricProcess.BuildPhysicsTable(*gamma);
  } 
  
  else if(processID == 2) {

    ComptonProcess.BuildPhysicsTable(*gamma);
  } 
    
  else if(processID == 3) {

    RayleighProcess.BuildPhysicsTable(*gamma);
  } 

  else if(processID == 4) {

    GammaConvProcess.BuildPhysicsTable(*gamma);
  } 
  else if(processID == 5) {

    cout<<"****** BR BuildPhysicsTable:Table already Built *******"<<endl;
  } 
  else if(processID == 6) {

    cout<<"****** IO BuildPhysicsTable:Table already Built *******"<<endl;
  } 

  else {
    
    G4Exception("No such process!\n");
  }

  // Mean FreePath Test
  G4Material* apttoMaterial ;
  G4String MaterialName ;
  ///***********************************************************************
  // UNCOMMENT LATER OR NEVER I HAVE IT IN 6 diff files
  //***********************************************************************
  G4double minArg = 100*eV,maxArg = 100*GeV, argStp;
  const G4int pntNum = 300;
  G4double Tkin[pntNum+1];
  G4double meanFreePath ;

  argStp = (log10(maxArg)-log10(minArg))/pntNum;
  
  for(G4int d = 0; d < pntNum+1; d++){ 
    
    Tkin[d] = pow(10,(log10(minArg) + d*argStp));
  }
  
  for ( G4int J = 0 ; J < theMaterialTable->length() ; J++ ){

    apttoMaterial = (*theMaterialTable)[ J ] ;
    MaterialName  = apttoMaterial->GetName() ;

    LogicalFrame->SetMaterial(apttoMaterial); 

    for (G4int i=0 ; i<pntNum; i++){

      if(processID != 5 && processID != 6)
	photon.SetKineticEnergy(Tkin[i]);
      else
	elecT.SetKineticEnergy(Tkin[i]);

      if(processID == 1){

	meanFreePath = PhotoElectricProcess.GetMeanFreePath(aTrack, 1., condition);
      }
      
      else if(processID == 2) {
	
	meanFreePath = ComptonProcess.GetMeanFreePath(aTrack, 1., condition);
      } 
      
      else if(processID == 3) {
	
	meanFreePath = RayleighProcess.GetMeanFreePath(aTrack, 1., condition);
      } 
      
      else if(processID == 4) {
	
	meanFreePath = GammaConvProcess.GetMeanFreePath(aTrack, 1., condition);
      } 
      
      else if(processID == 5) {
	
	meanFreePath = BremstrahlungProcess.GetMeanFreePath(aTrack, 1., condition);
      } 
      
      else if(processID == 6) {
	
	meanFreePath = IonisationProcess.GetMeanFreePath(aTrack, 1., condition);
      } 
      
      else {
	
	G4Exception("No such process!\n");
      }

       ntuple3->column("matind",J);
       ntuple3->column("kinen",Tkin[i]);
       ntuple3->column("mfp",meanFreePath);
       ntuple3->dumpData();

      //      outFile3<<setw(4)<<J<<setw(14)<<Tkin[i]<<setw(14)<<meanFreePath<<endl;    
    }
  }// for loop on materials
  //END OF COMMENT */
  // --------- Test the DoIt for the Photon Absorption

  apttoMaterial = (*theMaterialTable)(imat) ;
  
  LogicalFrame->SetMaterial(apttoMaterial); 
  if(processID != 5 && processID != 6){
    photon.SetKineticEnergy(InitEnergy*MeV);
    photon.SetMomentumDirection(InitX, InitY, InitZ);
  }
  else{
    elecT.SetKineticEnergy(InitEnergy*MeV);
    elecT.SetMomentumDirection(InitX, InitY, InitZ);
  }
  // PostStepDoIt calls
  G4int iteration = 0;   
 
  G4VParticleChange* adummy;
  G4Track* aFinalParticle;
  G4String aParticleName;  

  do {
    
    if(processID == 1){

     adummy = PhotoElectricProcess.PostStepDoIt(aTrack, aStep);

    } 

    else if(processID == 2) {

      adummy = ComptonProcess.PostStepDoIt(aTrack, aStep);
    } 
    
    else if(processID == 3) {

      adummy = RayleighProcess.PostStepDoIt(aTrack, aStep);
    } 
    
    else if(processID == 4) {

      adummy = GammaConvProcess.PostStepDoIt(aTrack, aStep);
    } 
    
    else if(processID == 5) {

      adummy = BremstrahlungProcess.PostStepDoIt(aTrack, aStep);
    } 
    
    else if(processID == 6) {

      adummy = IonisationProcess.PostStepDoIt(aTrack, aStep);
    } 
    
    else {

      G4Exception("No such process!\n");
    }

    G4ParticleChange* aParticleChange = (G4ParticleChange*) adummy;
  
    // ------------ book primary physical quantities -------------
    G4double pEnChange = 0, pxChange = 0, pyChange = 0, pzChange = 0, PChange = 0;

    pEnChange = aParticleChange->GetEnergyChange();
    pxChange  = aParticleChange->GetMomentumChange()->x();
    pyChange  = aParticleChange->GetMomentumChange()->y();
    pzChange  = aParticleChange->GetMomentumChange()->z();
    PChange   = sqrt(pxChange*pxChange+pyChange*pyChange+pzChange*pzChange);

    // ---- secondaries histos ----    
    G4cout<<"E and p of the primary particle: "<<pEnChange<<"  "<<pxChange<<"  "
	  <<pyChange<<"  "<<pzChange<<endl;

    ntuple1->column("ench", pEnChange);
    ntuple1->column("pxch", pxChange);
    ntuple1->column("pych", pyChange);
    ntuple1->column("pzch", pzChange);
    ntuple1->column("pch", PChange);
    
    //    outFile1<<setw(13)<<pEnChange<<setw(13)<<pxChange<<setw(13)
    //	    <<pyChange<<setw(13)<<pzChange<<endl;
    
    ntuple1->dumpData(); 
    
    // ------------ book secondaries physical quantities ---------
    
    G4double e = 0, eKin = 0, Px = 0, Py = 0, Pz = 0, P = 0, ShID = 0;
    
    //    hNSec->accumulate(aParticleChange->GetNumberOfSecondaries());
    // hDebug->accumulate(aParticleChange->GetLocalEnergyDeposit());
    
    for (G4int i = 0; i < (aParticleChange->GetNumberOfSecondaries()); i++) {

      // The following two items should be filled per event, not
      // per secondary; filled here just for convenience, to avoid
      // complicated logic to dump ntuple when there are no secondaries

       ntuple2->column("nsec",aParticleChange->GetNumberOfSecondaries());
       ntuple2->column("deposit",aParticleChange->GetLocalEnergyDeposit());
      
      aFinalParticle = aParticleChange->GetSecondary(i) ;
      
      e    = aFinalParticle->GetTotalEnergy();
      eKin = aFinalParticle->GetKineticEnergy();
      Px   = (aFinalParticle->GetMomentum()).x();
      Py   = (aFinalParticle->GetMomentum()).y();
      Pz   = (aFinalParticle->GetMomentum()).z();
      P    = sqrt(Px*Px+Py*Py+Pz*Pz);
      if(processID == 1){

	ShID = PhotoElectricProcess.GetTransitionShell(i);
      }
      else if(processID == 6){

	ShID = IonisationProcess.GetTransitionShell(i);
      }

      aParticleName = aFinalParticle->GetDefinition()->GetParticleName();
      G4cout<<aParticleName<<": "
	  <<" "<<e<<"  "<<eKin<<"  "<<Px<<"  "<<Py<<"  "<<Pz<<" ***"<<endl;      
      hEKin->accumulate(eKin);
      hP->accumulate(sqrt(Px*Px+Py*Py+Pz*Pz));

      G4int ptype;
      if(aParticleName == "gamma") ptype = 0;
      else if(aParticleName == "e-") ptype = -1;
      else if(aParticleName == "e+") ptype = 1;

      // Fill the secondaries ntuple
      ntuple2->column("px", Px);
      ntuple2->column("py", Py);
      ntuple2->column("pz", Pz);
      ntuple2->column("p", P);
      ntuple2->column("e", e);
      ntuple2->column("ekin", eKin);
      ntuple2->column("ptype", ptype);
      ntuple2->column("sh", ShID);

      ntuple2->dumpData(); 

      // Print secondaries on a file
      //      outFile2<<setw(3)<<aParticleChange->GetNumberOfSecondaries()<<setw(14)
      //	      <<aParticleChange->GetLocalEnergyDeposit()<<setw(14)<<e<<setw(14)
      //	      <<eKin<<setw(14)<<sqrt(Px*Px+Py*Py+Pz*Pz)<<setw(14)<<Px<<setw(14)
      //	      <<Py<<setw(14)<<Pz<<setw(3)<<ptype<<endl;

      delete aParticleChange->GetSecondary(i);
    }
    
    aParticleChange->Clear();
    iteration++;
    G4cout << "******* Iteration = " << iteration << endl;
    
  }  while (iteration < niter) ;
  cout<<"Iteration number: "<<endl;
  hbookManager->write();
  delete hbookManager;

  // delete materials and elements
  delete Be;
  delete Graphite;
  delete Al;
  delete LAr;
  delete Fe;
  delete Cu;
  delete W;
  delete Pb;
  delete U;
  delete H;
  delete maO;
  delete C;
  delete Cs;
  delete I;
  delete O;
  delete water;
  delete ethane;
  delete csi;
  delete Step;
  delete touche;
  //  outFile1.close();
  //  outFile2.close();
  //  outFile3.close();

  cout<<"END OF THE MAIN PROGRAM"<<endl;
}













