// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4hLETestAlongPostStep.cc,v 1.1 2000-05-03 16:32:07 chauvie Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      File name:     G4hLETestAlongPostStep.cc
//
//      Author:        Stephane Chauvie 
// 
//      Creation date: 2 May 2000
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>

#include "G4Material.hh"
#include "G4ProcessManager.hh"
#include "G4hLowEnergyIonisation.hh"
#include "G4EnergyLossTables.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"

#include "G4Step.hh"
#include "G4GRSVolume.hh"

#include "G4UnitsTable.hh"
#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"

HepTupleManager* hbookManager;

main()
{

  // Setup

  G4int niter=100000;
  G4int imat=3;
  G4int particleID=1;
  G4int test=0;

  G4cout << "Test AlongStepDoIt [1] or PostStepDoIt [2] ?" << G4endl;
  cin >> test;
  G4cout << "How many interactions? [100000], Which material? [3]" << G4endl;
  cin >> niter >> imat ;

  G4cout <<"Which particle?"<<endl<<setw(60)<<"[1] = Proton, [2] =AntiProton"<< G4endl;
  cin>>particleID;
  
  G4double InitEnergy = 1*MeV, InitX = 0., InitY = 0., InitZ = 1.;
  G4cout<<"Enter the initial particle energy E"<<G4endl; 
  G4cin >> InitEnergy ;

  G4cout.setf( ios::scientific, ios::floatfield );
  // -------------------------------------------------------------------

  // ---- HBOOK initialization


  hbookManager = new HBookFile("G4hLETestAlongPostStep.hbook", 58);
  assert (hbookManager != 0);
  
  // ---- Book a histogram and ntuples
  G4cout<<"Hbook file name: "<<((HBookFile*) hbookManager)->filename()<<endl;
  
  // ---- primary ntuple ------
  HepTuple* ntuple1 = hbookManager->ntuple("Primary Ntuple");
  assert (ntuple1 != 0);
  
  // ---- secondary ntuple ------
  HepTuple* ntuple2 = hbookManager->ntuple("Secondary Ntuple");
  assert (ntuple2 != 0);
   
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
  G4Material* Si  = new G4Material("Silicon",   14., 28.055*g/mole, 2.33*g/cm3);
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

  G4double dimx = 1*mm, dimy = 1*mm, dimz = 1*mm;
  
  // Geometry definitions
  G4Box* theFrame = new G4Box ("Frame",dimx, dimy, dimz);
  
  G4LogicalVolume* LogicalFrame = new G4LogicalVolume(theFrame,
						      (*theMaterialTable)(imat),
						      "LFrame", 0, 0, 0);
  
  G4PVPlacement* PhysicalFrame = new G4PVPlacement(0,G4ThreeVector(),
						   "PFrame",LogicalFrame,0,false,0);
  
  // the center-of-mass of the cube should be located at the origin!

  //--------- Particle definition ---------

  G4ParticleDefinition* gamma = G4Gamma::GammaDefinition();
  G4ParticleDefinition* electron = G4Electron::ElectronDefinition();
  G4ParticleDefinition* proton = G4Proton::ProtonDefinition();
  G4ParticleDefinition* antiproton = G4AntiProton::AntiProtonDefinition();
  
  electron->SetCuts(1e-3*mm);
  proton->SetCuts(1e-3*mm);
  antiproton->SetCuts(1e-3*mm);
  
  //--------- Processes definition ---------

  G4hLowEnergyIonisation* hIonisationProcess = new G4hLowEnergyIonisation;
  	
	G4ProcessManager* theProtonProcessManager = new G4ProcessManager(proton);
  	proton->SetProcessManager(theProtonProcessManager);
	theProtonProcessManager->AddProcess(hIonisationProcess);
  	//theProtonProcessManager->SetProcessOrdering(hIonisationProcess,idxAlongStep,1);
  	//theProtonProcessManager->SetProcessOrdering(hIonisationProcess,idxPostStep,1);

  	G4ProcessManager* theAntiProtonProcessManager = new G4ProcessManager(antiproton);
  	antiproton->SetProcessManager(theAntiProtonProcessManager);
        theAntiProtonProcessManager->AddProcess(hIonisationProcess);
  	//theAntiProtonProcessManager->SetProcessOrdering(hIonisationProcess,idxAlongStep,1);
  	//theAntiProtonProcessManager->SetProcessOrdering(hIonisationProcess,idxPostStep,1);
        hIonisationProcess->SetAntiProtonStoppingOn();
  
  // -------- create 1 Dynamic Particle  ----

  G4double pEnergy = InitEnergy*MeV;

  G4ParticleMomentum pDirection(InitX,InitY,InitZ);
  
  G4DynamicParticle p(G4Proton::Proton(),pDirection,pEnergy);
  G4DynamicParticle pbar(G4AntiProton::AntiProton(),pDirection,pEnergy);
  

  //--------- track definition (for this test ONLY!)------------

  G4ThreeVector aPosition(0.,0.,0.);
  G4ThreeVector aPositionf(0.,0.,0.001*mm);
  G4double aTime = 0. ;

  G4Track* ptrack;
  G4Track* ptrackbar;

  ptrack = new G4Track(&p,aTime,aPosition) ;
  ptrackbar = new G4Track(&pbar,aTime,aPosition) ;

  G4Track& aTrack = (*ptrack);
  if(particleID==2) aTrack = (*ptrackbar) ;
  
  // do I really need this?

  G4GRSVolume* touche = new G4GRSVolume(PhysicalFrame, NULL, aPosition);   
  ptrack->SetTouchable(touche);
  ptrackbar->SetTouchable(touche);
 
 // -------- create 1 Step (for this test only)----  

  G4Step* Step = new G4Step();  
  G4Step& aStep = (*Step);
  Step->SetTrack(ptrack);
  if(particleID==2)Step->SetTrack(ptrackbar);
  
  G4StepPoint* aPoint = new G4StepPoint();
  (*aPoint).SetPosition(aPosition);
  G4double safety = 10000.*cm;
  (*aPoint).SetSafety(safety);
  (*Step).SetPreStepPoint(aPoint);
  //(*aPoint).SetPosition(aPositionf);
  //(*Step).SetPostStepPoint(aPoint);
  
  // --------- check applicability
  
  G4ParticleDefinition* ProtonDefinition = p.GetDefinition();
  G4ParticleDefinition* AntiProtonDefinition = pbar.GetDefinition();

  if(! (hIonisationProcess->IsApplicable(*ProtonDefinition)
	|| hIonisationProcess->IsApplicable(*AntiProtonDefinition)) )
    {
      G4Exception("FAIL: *** Not Applicable ***\n");
    }

  // Initialize the physics tables for ALL processes

if( particleID == 1){
  	
	hIonisationProcess->BuildPhysicsTable(*ProtonDefinition);
	
  } else {
  	
	hIonisationProcess->BuildPhysicsTable(*AntiProtonDefinition);
	
  }	
  
  G4Material* apttoMaterial ;
  G4String MaterialName ;


  // --------- Test the DoIt for the hIonization

  apttoMaterial = (*theMaterialTable)(imat) ;
  
  LogicalFrame->SetMaterial(apttoMaterial); 

  G4int iteration = 0;   
 
  G4VParticleChange* adummy;
  G4Track* aFinalParticle;
  G4String aParticleName;  
  G4double dedx=0,dedxnow=0,dx=0;
  G4ThreeVector vecdx;
  G4double xChange, yChange, zChange;
  do {
    (*Step).SetStepLength(1*micrometer);
    G4cout<<"Step Lenght ="<<(*Step).GetStepLength()/mm<<endl;
    aTrack.SetStep(Step); //this function should be added because
                          //the Step Length is not accessible from
			  //aTrack but you should use Step
    G4cout<<aTrack.GetStep()->GetStepLength()/mm<<endl;
    if(particleID==1) hIonisationProcess->GetConstraints(&p,apttoMaterial);
    if(particleID==2) hIonisationProcess->GetConstraints(&pbar,apttoMaterial);
    if(test==1)adummy = hIonisationProcess->AlongStepDoIt(aTrack, aStep);
    if(test==2)adummy = hIonisationProcess->PostStepDoIt(aTrack, aStep);
    G4ParticleChange* aParticleChange = (G4ParticleChange*) adummy;
  
    // ------------ book primary physical quantities -------------
    G4double pEnChange = 0, pxChange = 0, pyChange = 0, pzChange = 0, pChange = 0;

    pEnChange = aParticleChange->GetEnergyChange();
    dedx =InitEnergy-pEnChange ;
    dedxnow= dedx/(*Step).GetStepLength();
    
    pxChange  = aParticleChange->GetMomentumChange()->x();
    pyChange  = aParticleChange->GetMomentumChange()->y();
    pzChange  = aParticleChange->GetMomentumChange()->z();
    pChange   = sqrt(pxChange*pxChange+pyChange*pyChange+pzChange*pzChange);
    
    xChange=aParticleChange->GetPositionChange()->x();
    yChange=aParticleChange->GetPositionChange()->y();
    zChange=aParticleChange->GetPositionChange()->z();
    G4cout<<"Particle Position (x,y,z) = "<<xChange<<"  "<<yChange<<"   "<<zChange<<"   "<<endl;
    G4cout<<"Particle's Energy after the Step ="<<pEnChange<<endl;
    G4cout<<"Energy Loss (dE) ="<<dedx<<endl;
    G4cout<<"Stopping Power (dE/dx)="<<dedxnow<<endl;
    
    // ---- secondaries histos ----    
    G4cout<<"E and p of the primary particle: "<<pEnChange<<"  "<<pxChange<<"  "
	  <<pyChange<<"  "<<pzChange<<endl;

    ntuple1->column("energyf", pEnChange);
    ntuple1->column("de", dedx);
    ntuple1->column("dedx", dedxnow);
    ntuple1->column("pxch", pxChange);
    ntuple1->column("pych", pyChange);
    ntuple1->column("pzch", pzChange);
    ntuple1->column("pch",  zChange);  //prima era dx
    ntuple1->dumpData(); 
    
    // ------------ book secondaries physical quantities ---------
    
    G4double e = 0, eKin = 0, Px = 0, Py = 0, Pz = 0, P = 0, ShID = 0;
    
    hNSec->accumulate(aParticleChange->GetNumberOfSecondaries());
    hDebug->accumulate(aParticleChange->GetLocalEnergyDeposit());
    
    for (G4int i = 0; i < (aParticleChange->GetNumberOfSecondaries()); i++) {

      // The following two items should be filled per event, not
      // per secondary; filled here just for convenience, to avoid
      // complicated logic to dump ntuple when there are no secondaries

      aFinalParticle = aParticleChange->GetSecondary(i) ;
      
      e    = aFinalParticle->GetTotalEnergy();
      eKin = aFinalParticle->GetKineticEnergy();
      Px   = (aFinalParticle->GetMomentum()).x();
      Py   = (aFinalParticle->GetMomentum()).y();
      Pz   = (aFinalParticle->GetMomentum()).z();
      P    = sqrt(Px*Px+Py*Py+Pz*Pz);

      aParticleName = aFinalParticle->GetDefinition()->GetParticleName();
      G4cout << aParticleName << ": "
	     << " " << e << "  " 
	     << eKin << "  " 
	     << Px << "  " 
	     << Py << "  "
	     << Pz << " ***" << endl;   
   
      hEKin->accumulate(eKin);
      hP->accumulate(sqrt(Px*Px+Py*Py+Pz*Pz));

      G4int ptype;
      if(aParticleName == "proton") ptype = 0;
      else if(aParticleName == "e-") ptype = -1;
      else if(aParticleName == "e+") ptype = 1;
      else if(aParticleName == "antiproton") ptype = 2;

      // Fill the secondaries ntuple
      ntuple2->column("px", Px);
      ntuple2->column("py", Py);
      ntuple2->column("pz", Pz);
      ntuple2->column("p", P);
      ntuple2->column("e", e);
      ntuple2->column("ekin", eKin);
      ntuple2->column("ptype", ptype);

      ntuple2->dumpData(); 

      delete aParticleChange->GetSecondary(i);
    }
    
    aParticleChange->Clear();
    iteration++;
    G4cout << "******* Iteration = " << iteration << endl;
    
  }  while (iteration < niter) ;

  cout <<"Iteration number: " << endl;
  hbookManager->write();
  delete hbookManager;

  // delete materials and elements
  delete Be;
  delete Graphite;
  delete Al;
  delete Si;
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

  cout<<"END OF THE MAIN PROGRAM"<<endl;
}












