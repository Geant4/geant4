// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyGammaConversionTest.cc,v 1.1 2001-01-08 16:22:48 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      File name:     G4LowEnergyGammaConversionTest.cc
//
//      Author:        Francesco Longo 
// 
//      Creation date: 04 january 2001
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4VDiscreteProcess.hh"
#include "G4VProcess.hh"
#include "G4ProcessManager.hh"

#include "G4LowEnergyGammaConversion.hh"
#include "G4GammaConversion.hh"

#include "G4EnergyLossTables.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"

#include "G4LowEnergyBremsstrahlung.hh"
#include "G4LowEnergyIonisation.hh"
#include "G4eIonisation.hh"
#include "G4MultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

//#include "G4ComptonScattering.hh"
//#include "G4PhotoElectricEffect.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

#include "G4GRSVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"

#include "G4UnitsTable.hh"

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"


HepTupleManager* hbookManager;

G4int main()
{

  // Setup

  G4int nIterations = 100000;
  G4int materialId = 3;

  G4cout.setf( ios::scientific, ios::floatfield );

  // -------------------------------------------------------------------

  // ---- HBOOK initialization


  hbookManager = new HBookFile("gammatest.hbook", 58);
  assert (hbookManager != 0);
  
  // ---- Book a histogram and ntuples
  G4cout<<"Hbook file name: "<<((HBookFile*) hbookManager)->filename()<<endl;
  
  // ---- primary ntuple ------
  HepTuple* ntuple1 = hbookManager->ntuple("Primary Ntuple");
  assert (ntuple1 != 0);
  
  // ---- secondary ntuple ------
  HepTuple* ntuple2 = hbookManager->ntuple("Secondary Ntuple");
  assert (ntuple2 != 0);

  // ---- table ntuple ------
  HepTuple* ntuple3 = hbookManager->ntuple("Mean Free Path Ntuple");
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

  G4Material* Si  = new G4Material("Silicon",   14., 28.055*g/mole, 2.33*g/cm3);
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


  // Interactive set-up

  G4cout << "How many interactions? " << G4endl;
  G4cin >> nIterations;

  if (nIterations <= 0) G4Exception("Wrong input");

  G4double initEnergy = 1*MeV; 
  G4double initX = 0.; 
  G4double initY = 0.; 
  G4double initZ = 1.;
  
  G4cout << "Enter the initial particle energy E (MeV)" << G4endl; 
  G4cin >> initEnergy ;
  
  initEnergy = initEnergy * MeV;
  
  if (initEnergy  <= 0.) G4Exception("Wrong input");

  static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  G4int nMaterials = theMaterialTable->length();

  G4cout << "Available materials are: " << G4endl;
  for (G4int mat = 0; mat < nMaterials; mat++)
    {
      G4cout << mat << ") "
	     << (*theMaterialTable)(mat)->GetName()
	     << G4endl;
    }
  
  G4cout << "Which material? " << G4endl;
  G4cin >> materialId;
  
  G4Material* material = (*theMaterialTable)(materialId) ;

  G4cout << "The selected material is: "
	 << material->GetName()
	 << G4endl;
  
  G4double dimX = 1*mm;
  G4double dimY = 1*mm;
  G4double dimZ = 1*mm;
  
  // Geometry 
  
  G4Box* theFrame = new G4Box ("Frame",dimX, dimY, dimZ);
  
  G4LogicalVolume* logicalFrame = new G4LogicalVolume(theFrame,
						      (*theMaterialTable)(materialId),
						      "LFrame", 0, 0, 0);
  logicalFrame->SetMaterial(material); 
  
  G4PVPlacement* physicalFrame = new G4PVPlacement(0,G4ThreeVector(),
						   "PFrame",logicalFrame,0,false,0);
  
  // Particle definitions
  
  G4ParticleDefinition* gamma = G4Gamma::GammaDefinition();
  G4ParticleDefinition* electron = G4Electron::ElectronDefinition();
  G4ParticleDefinition* positron = G4Positron::PositronDefinition();
  
  gamma->SetCuts(1*micrometer);
  electron->SetCuts(1*micrometer);
  positron->SetCuts(1*micrometer);

  G4Gamma::SetEnergyRange(2.5e-4*MeV,1e5*MeV);
  G4Electron::SetEnergyRange(2.5e-4*MeV,1e5*MeV);
  G4Positron::SetEnergyRange(2.5e-4*MeV,1e5*MeV);

  G4cout<<"the cut in energy for gamma in: "<<
    (*theMaterialTable)(materialId)->GetName()
	<<" is: "<<G4Gamma::GetCutsInEnergy()[materialId]<<G4endl;
  G4cout<<"the cut in energy for e- in: "<<
    (*theMaterialTable)(materialId)->GetName()
	<<" is: "<<G4Electron::GetCutsInEnergy()[materialId]<<G4endl;
  
  // Processes 
  
  
  G4int processType;
  G4cout << "LowEnergy [1] or Standard [2] Gamma Conversion?" << G4endl;
  G4cin >> processType;
  if ( !(processType == 1 || processType == 2))
    {
      G4Exception("Wrong input");
    }

  G4VDiscreteProcess* gammaProcess;

  if (processType == 1)
    {
      gammaProcess = new G4LowEnergyGammaConversion();
    }
  else
    {
      gammaProcess = new G4GammaConversion();
    }
  
  
  G4VProcess* theeminusMultipleScattering = new G4MultipleScattering();
  G4VProcess* theeminusIonisation         = new G4eIonisation();
  G4VProcess* theeminusBremsstrahlung     = new G4eBremsstrahlung();
  G4VProcess* theeplusMultipleScattering  = new G4MultipleScattering();
  G4VProcess* theeplusIonisation          = new G4eIonisation();
  G4VProcess* theeplusBremsstrahlung      = new G4eBremsstrahlung();
  G4VProcess* theeplusAnnihilation        = new G4eplusAnnihilation();

  //----------------
  // process manager  
  //----------------

  // gamma
  
  G4ProcessManager* gProcessManager = new G4ProcessManager(gamma);
  gamma->SetProcessManager(gProcessManager);
  gProcessManager->AddDiscreteProcess(gammaProcess);
  G4ForceCondition* condition;

  //electron
  
  G4ProcessManager* eProcessManager = new G4ProcessManager(electron);
  electron->SetProcessManager(eProcessManager);
  eProcessManager->AddProcess(theeminusMultipleScattering);
  eProcessManager->AddProcess(theeminusIonisation);
  eProcessManager->AddProcess(theeminusBremsstrahlung);
  
  //positron
  
  G4ProcessManager* pProcessManager = new G4ProcessManager(positron);
  positron->SetProcessManager(pProcessManager);
  pProcessManager->AddProcess(theeplusMultipleScattering);
  pProcessManager->AddProcess(theeplusIonisation);
  pProcessManager->AddProcess(theeplusBremsstrahlung);
  pProcessManager->AddProcess(theeplusAnnihilation);
  
  //--------------
  // set ordering   
  //--------------


  eProcessManager->
    SetProcessOrdering(theeminusMultipleScattering, idxAlongStep,1);
  eProcessManager->
    SetProcessOrdering(theeminusIonisation,         idxAlongStep,2);
      
  eProcessManager->
    SetProcessOrdering(theeminusMultipleScattering, idxPostStep,1);
  eProcessManager->
    SetProcessOrdering(theeminusIonisation,         idxPostStep,2);
  eProcessManager->
    SetProcessOrdering(theeminusBremsstrahlung,     idxPostStep,3);



  pProcessManager->SetProcessOrderingToFirst(theeplusAnnihilation, idxAtRest);
  pProcessManager->
    SetProcessOrdering(theeplusMultipleScattering, idxAlongStep,1);
  pProcessManager->
    SetProcessOrdering(theeplusIonisation,         idxAlongStep,2);

  pProcessManager->
    SetProcessOrdering(theeplusMultipleScattering, idxPostStep,1);
  pProcessManager->
    SetProcessOrdering(theeplusIonisation,         idxPostStep,2);
  pProcessManager->
    SetProcessOrdering(theeplusBremsstrahlung,     idxPostStep,3);
  pProcessManager->
    SetProcessOrdering(theeplusAnnihilation,       idxPostStep,4);
  
  // G4LowEnergyIonisation IonisationProcess;
  // eProcessManager->AddProcess(&IonisationProcess);
  // eProcessManager->SetProcessOrdering(&IonisationProcess,idxAlongStep,1);
  // eProcessManager->SetProcessOrdering(&IonisationProcess,idxPostStep, 1);
  
  // G4LowEnergyBremsstrahlung BremstrahlungProcess;
  // eProcessManager->AddProcess(&BremstrahlungProcess);
  // eProcessManager->SetProcessOrdering(&BremstrahlungProcess,idxAlongStep,1);
  // eProcessManager->SetProcessOrdering(&BremstrahlungProcess,idxPostStep, 1);
  
  // G4eIonisation IonisationPlusProcess;
  // pPositronProcessManager->AddProcess(&IonisationPlusProcess);
  // pProcessManager->
  //        SetProcessOrdering(&IonisationPlusProcess,idxAlongStep,1);
  // pProcessManager->SetProcessOrdering(&IonisationPlusProcess,idxPostStep,1);



  // Create a DynamicParticle  
  
  G4double eEnergy = initEnergy*MeV;
  G4ParticleMomentum eDirection(initX,initY,initZ);
  G4DynamicParticle dynamicGamma(G4Gamma::Gamma(),eDirection,eEnergy);

  dynamicGamma.DumpInfo(0);
  
  // Track 

  G4ThreeVector aPosition(0.,0.,0.);
  G4double aTime = 0. ;
  
  G4Track* gTrack = new G4Track(&dynamicGamma,aTime,aPosition);

  G4GRSVolume* touche = new G4GRSVolume(physicalFrame, NULL, aPosition);   
  gTrack->SetTouchable(touche);


  // Step 

  G4Step* step = new G4Step();  
  step->SetTrack(gTrack);

  G4StepPoint* aPoint = new G4StepPoint();
  aPoint->SetPosition(aPosition);
  aPoint->SetMaterial(material);
  G4double safety = 10000.*cm;
  aPoint->SetSafety(safety);
  step->SetPreStepPoint(aPoint);
  
  // Check applicability
  
  if (! (gammaProcess->IsApplicable(*gamma)))
    {
      G4Exception("Not Applicable");
    }
  else 
    {
      G4cout<< "applicability OK" << endl;
    }
  
  // Initialize the physics tables (in which material?)

  gammaProcess->BuildPhysicsTable(*gamma);

  theeminusMultipleScattering->BuildPhysicsTable(*electron);
  theeminusIonisation->BuildPhysicsTable(*electron);        
  theeminusBremsstrahlung->BuildPhysicsTable(*electron);
  theeplusMultipleScattering->BuildPhysicsTable(*positron);
  theeplusIonisation->BuildPhysicsTable(*positron);
  theeplusBremsstrahlung->BuildPhysicsTable(*positron);     
  theeplusAnnihilation->BuildPhysicsTable(*positron) ;

  G4cout<< "table OK" << endl;
  
  // Test GetMeanFreePath()
  
  G4Material* apttoMaterial ;
  G4String MaterialName ;
  
  G4double minArg = 100*eV,maxArg = 100*GeV, argStp;
  const G4int pntNum = 300;
  G4double Tkin[pntNum+1];
  G4double meanFreePath=0. ;

  argStp = (log10(maxArg)-log10(minArg))/pntNum;
  
  for(G4int d = 0; d < pntNum+1; d++)
    { 
      Tkin[d] = pow(10,(log10(minArg) + d*argStp));
    }
 
  G4double sti = 1.*mm;
  step->SetStepLength(sti);
  
  //  for ( G4int J = 0 ; J < nMaterials ; J++ )
  //  {
  apttoMaterial = (*theMaterialTable)(materialId) ;
  MaterialName  = apttoMaterial->GetName() ;
  logicalFrame->SetMaterial(apttoMaterial); 
  
  gTrack->SetStep(step);

  G4LowEnergyGammaConversion* gammaLowEProcess =
    (G4LowEnergyGammaConversion*) gammaProcess;
  G4GammaConversion* gammaStdProcess =
    (G4GammaConversion*) gammaProcess;
  
  
  for (G4int i=0 ; i<pntNum; i++)
    {
      dynamicGamma.SetKineticEnergy(Tkin[i]);
      if (processType == 1)
	{
	  meanFreePath=gammaLowEProcess
	    ->GetMeanFreePath(*gTrack, sti, condition);
	}
      else
	{
	  meanFreePath=gammaStdProcess
	    ->GetMeanFreePath(*gTrack, sti, condition);
	}

      ntuple3->column("kinen",Tkin[i]);
      ntuple3->column("mfp",meanFreePath/cm);
      ntuple3->dumpData();
    
      //      G4cout << meanFreePath/cm << G4endl;

    }
  G4cout << "Mean Free Path OK" << G4endl;
  
  // --------- Test the DoIt 
  
  G4cout << "DoIt in " << material->GetName() << G4endl;


  dynamicGamma.SetKineticEnergy(eEnergy);
  for (G4int iter=0; iter<nIterations; iter++)
    {
      
      step->SetStepLength(1*micrometer);
      
      G4cout  <<  "Iteration = "  <<  iter 
	      << "  -  Step Length = " 
	      << step->GetStepLength()/mm << " mm "
	      << G4endl;
      
      gTrack->SetStep(step); 
 
      //      G4cout  <<  "Iteration = "  <<  iter 
      //	      << "  -  Step Length = " 
      //      << step->GetStepLength()/mm << " mm "
      //      << G4endl;
      
      //G4cout << gTrack->GetStep()->GetStepLength()/mm 
      //     << G4endl;
      
      G4VParticleChange* dummy;
      dummy = gammaProcess->PostStepDoIt(*gTrack, *step);
      
      G4ParticleChange* particleChange = (G4ParticleChange*) dummy;
      
      // Primary physical quantities 
      
      G4double energyChange = particleChange->GetEnergyChange();
      
      G4double dedx = initEnergy - energyChange ;
      G4double dedxNow = dedx / (step->GetStepLength());
      
      G4ThreeVector eChange = 
	particleChange->CalcMomentum(energyChange,
				     (*particleChange->GetMomentumChange()),
				     particleChange->GetMassChange());
      
      G4double pxChange  = eChange.x();
      G4double pyChange  = eChange.y();
      G4double pzChange  = eChange.z();
      G4double pChange   = 
	sqrt(pxChange*pxChange + pyChange*pyChange + pzChange*pzChange);
      
      G4double xChange = particleChange->GetPositionChange()->x();
      G4double yChange = particleChange->GetPositionChange()->y();
      G4double zChange = particleChange->GetPositionChange()->z();
      G4double thetaChange = particleChange->GetPositionChange()->theta();
      
      G4cout << "---- Primary after the step ---- " << G4endl;
 
      //      G4cout << "Position (x,y,z) = " 
      //	     << xChange << "  " 
      //	     << yChange << "   " 
      //	     << zChange << "   " 
      //	     << G4endl;

      G4cout << "---- Energy: " << energyChange/MeV << " MeV,  " 
	     << "(px,py,pz): ("
	     << pxChange/MeV << ","
	     << pyChange/MeV << "," 
	     << pzChange/MeV << ") MeV"
	     << G4endl;
      
      G4cout << "---- Energy loss (dE) = " << dedx/keV << " keV" << G4endl;
      //      G4cout << "Stopping power (dE/dx)=" << dedxNow << G4endl;
      
      // Secondaries 

      ntuple1->column("eprimary", initEnergy);
      ntuple1->column("energyf", energyChange);
      ntuple1->column("de", dedx);
      ntuple1->column("dedx", dedxNow);
      ntuple1->column("pxch", xChange);
      ntuple1->column("pych", pyChange);
      ntuple1->column("pzch", pzChange);
      ntuple1->column("pch", zChange);  
      ntuple1->column("thetach", thetaChange);  
      ntuple1->dumpData(); 

      // Secondaries physical quantities 
      
      hNSec->accumulate(particleChange->GetNumberOfSecondaries());
      hDebug->accumulate(particleChange->GetLocalEnergyDeposit());
      
      G4cout << " secondaries " << 
	particleChange->GetNumberOfSecondaries() << G4endl;
      
      for (G4int i = 0; i < (particleChange->GetNumberOfSecondaries()); i++) 
	{
	  // The following two items should be filled per event, not
	  // per secondary; filled here just for convenience, to avoid
	  // complicated logic to dump ntuple when there are no secondaries
	  
	  G4Track* finalParticle = particleChange->GetSecondary(i) ;
	  
	  G4double e    = finalParticle->GetTotalEnergy();
	  G4double eKin = finalParticle->GetKineticEnergy();
	  G4double px   = (finalParticle->GetMomentum()).x();
	  G4double py   = (finalParticle->GetMomentum()).y();
	  G4double pz   = (finalParticle->GetMomentum()).z();
	  G4double theta   = (finalParticle->GetMomentum()).theta();
	  G4double p   = sqrt(px*px+py*py+pz*pz);
	  
	  if (e > initEnergy)
	    {
	      G4cout << "WARNING: eFinal > eInit " << G4endl;
	      //	     << e
	      //		     << " > " initEnergy 
	      
	    }
	  
	  G4String particleName = 
	    finalParticle->GetDefinition()->GetParticleName();
	  G4cout  << "==== Final " 
		  <<  particleName  <<  " "  
		  << "energy: " <<  e/MeV  <<  " MeV,  " 
		  << "eKin: " <<  eKin/MeV  <<  " MeV, " 
		  << "(px,py,pz): ("
		  <<  px/MeV  <<  "," 
		  <<  py/MeV  <<  ","
		  <<  pz/MeV  << ") MeV "
		  <<  G4endl;   
	  
	  hEKin->accumulate(eKin);
	  hP->accumulate(p);
	  
	  G4int partType;
	  if (particleName == "e-") partType = 1;
	  else if (particleName == "e+") partType = 2;
	  else if (particleName == "gamma") partType = 3;
	  
	  // Fill the secondaries ntuple

          ntuple2->column("event",iter);
          ntuple2->column("eprimary",initEnergy);
	  ntuple2->column("px", px);
	  ntuple2->column("py", py);
	  ntuple2->column("pz", pz);
	  ntuple2->column("p", p);
	  ntuple2->column("e", e);
	  ntuple2->column("theta", theta);
	  ntuple2->column("ekin", eKin);
	  ntuple2->column("type", partType);
	  
	  ntuple2->dumpData(); 
	  
	  delete particleChange->GetSecondary(i);
	}
      
      particleChange->Clear();
      
    } 
  
  
  G4cout  << "Iteration number: "  <<  iter << G4endl;
  hbookManager->write();
  delete hbookManager;
  
  delete step;


  cout << "END OF THE MAIN PROGRAM" << G4endl;
}

















