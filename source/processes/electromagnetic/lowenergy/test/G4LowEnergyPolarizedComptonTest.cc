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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4LowEnergyPolarizedComptonTest.cc,v 1.4 2001-09-10 18:07:55 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4ComptonScatteringTest.cc
//
//      Author:        Francesco Longo & Gerardo Depaola
// 
//      Creation date: 23 january 2001
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

#include "G4ComptonScattering.hh"
#include "G4PolarizedComptonScattering.hh"
#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyPolarizedCompton.hh"

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

  //  G4cout.setf( ios::scientific, ios::floatfield );

  // -------------------------------------------------------------------

  // ---- HBOOK initialization


  hbookManager = new HBookFile("comptontest.hbook", 58);
  assert (hbookManager != 0);

  // ---- Book a histogram and ntuples

  G4cout<<"Hbook file name: "<<((HBookFile*) hbookManager)->filename()<<G4endl;
 

  G4double initEnergy = 1*MeV; 
  G4double initX = 0.; 
  G4double initY = 0.; 
  G4double initZ = 1.;
  
  G4cout << "Enter the initial particle energy E (keV)" << G4endl; 
  G4cin >> initEnergy ;
  initEnergy = initEnergy * keV;
  G4double limit = initEnergy/keV;

  G4cout << limit << G4endl;
  
  if (initEnergy  <= 0.) G4Exception("Wrong input");
 
  // ---- primary ntuple ------
  HepTuple* ntuple1 = hbookManager->ntuple("Primary Ntuple");
  assert (ntuple1 != 0);
  
  // ---- secondary ntuple ------
  HepTuple* ntuple2 = hbookManager->ntuple("Secondary Ntuple");
  assert (ntuple2 != 0);

  /*  
      // ---- table ntuple ------
      HepTuple* ntuple3 = hbookManager->ntuple("Mean Free Path Ntuple");
      assert (ntuple3 != 0);
  */
  
  // ---- secondaries histos ----

  HepHistogram* heETot;
  heETot = hbookManager->histogram("Electron Total Energy", 100,0.,limit);
  assert (heETot != 0);  
  
  HepHistogram* heP;
  heP = hbookManager->histogram("Electron Momentum", 100,0.,limit);
  assert (heP != 0);  

  HepHistogram* hgETot;
  hgETot = hbookManager->histogram("Gamma Total Energy", 100,0.,limit);
  assert (hgETot != 0);  
  
  HepHistogram* hgP;
  hgP = hbookManager->histogram("Gamma Momentum", 100,0.,limit);
  assert (hgP != 0);  

  HepHistogram* hgTheta;
  hgTheta = hbookManager->histogram("Theta Scattered Gamma ", 100,0.,4.);
  assert (hgTheta != 0);  

  HepHistogram* hgPhi;
  hgPhi = hbookManager->histogram("Phi Scattered Gamma ", 100,-4.,4.);
  assert (hgPhi != 0);  

  HepHistogram* hSumE;
  hSumE = hbookManager->histogram("Total Energy", 100,0.,2*limit);
  assert (hSumE != 0);  

  HepHistogram* hgRapp;
  hgRapp = hbookManager->histogram("Energy Theta Relation", 100,0.,2.);
  assert (hgRapp != 0);  

  HepHistogram* hNSec;
  hNSec = hbookManager->histogram("Number of secondaries", 100,0.,10.);
  assert (hNSec != 0);  

  HepHistogram* hDebug;
  hDebug = hbookManager->histogram("Debug", 100,0.,limit);
  assert (hDebug != 0);  
  

  //--------- Materials definition ---------

  G4Material* Si  = new G4Material("Silicon",   14., 28.055*g/mole, 2.33*g/cm3);
  G4Material* Fe  = new G4Material("Iron",      26., 55.85*g/mole, 7.87*g/cm3);
  G4Material* Cu  = new G4Material("Copper",    29., 63.55*g/mole, 8.96*g/cm3);
  G4Material*  W  = new G4Material("Tungsten", 74., 183.85*g/mole, 19.30*g/cm3);
  G4Material* Pb  = new G4Material("Lead",      82., 207.19*g/mole, 11.35*g/cm3);
  G4Element*   H  = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   O  = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
  G4Element*   C  = new G4Element ("Carbon"  , "C", 6. , 12.00*g/mole);
  G4Element*  Cs  = new G4Element ("Cesium"  , "Cs", 55. , 132.905*g/mole);
  G4Element*   I  = new G4Element ("Iodide"  , "I", 53. , 126.9044*g/mole);
  G4Material*  maO = new G4Material("Oxygen", 8., 16.00*g/mole, 1.1*g/cm3);

  G4Material* csi = new G4Material ("CsI" , 4.53*g/cm3, 2);
  csi->AddElement(Cs,1);
  csi->AddElement(I,1);


  // Interactive set-up

  G4cout << "How many interactions? " << G4endl;
  G4cin >> nIterations;

  if (nIterations <= 0) G4Exception("Wrong input");

 

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
  
  gamma->SetCuts(0.1*micrometer);
  electron->SetCuts(0.1*micrometer);
  positron->SetCuts(0.1*micrometer);

  G4Gamma::SetEnergyRange(2.5e-4*MeV,1e5*MeV);
  G4Electron::SetEnergyRange(2.5e-4*MeV,1e5*MeV);
  G4Positron::SetEnergyRange(2.5e-4*MeV,1e5*MeV);

  G4cout<<"the cut in energy for gamma in: "<<
    (*theMaterialTable)(materialId)->GetName()
	<<" is: "<<(G4Gamma::GetCutsInEnergy()[materialId])/keV
	<<" keV" <<G4endl;
  G4cout<<"the cut in energy for e- in: "<<
    (*theMaterialTable)(materialId)->GetName()
	<<" is: "<<(G4Electron::GetCutsInEnergy()[materialId])/keV
	<<" keV" <<G4endl;
  
  // Processes 
  

  G4int processType;
  G4cout << "LowEnergy [1] or Standard [2] Compton or Standard PolarizedCompton[3] or LowEnergyPolarizedCompton [4]" << G4endl;
  G4cin >> processType;
  if ( !(processType == 1 || processType == 2 
	 || processType == 3  || processType == 4))
    {
      G4Exception("Wrong input");
    }
  
  G4VDiscreteProcess* gammaProcess;
  

  if (processType == 1)
    {
      gammaProcess = new G4LowEnergyCompton();
    }
  else if (processType == 2)
    {
      gammaProcess = new G4ComptonScattering();
    }
  else if (processType == 3)
    {
      gammaProcess = new G4PolarizedComptonScattering();
    }
  else if (processType == 4)
    {
      gammaProcess = new G4LowEnergyPolarizedCompton();
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
  
  // Create a DynamicParticle  
  
  //  G4double eEnergy = initEnergy*keV;
  G4double eEnergy = initEnergy;
  
  //  G4cout << eEnergy/keV << " INIT ENERGY (keV)" << G4endl;

  G4ParticleMomentum eDirection(initX,initY,initZ);
  G4DynamicParticle dynamicGamma(G4Gamma::Gamma(),eDirection,eEnergy);

  G4cout << eDirection << " Direction" << G4endl;


  //  if (processType == 3 || processType == 4)
  //  {
      G4double PolX, PolY, PolZ;
      G4cout << "Polarization Vector" << G4endl;
      G4cin >> PolX >> PolY >> PolZ;
      dynamicGamma.SetPolarization(PolX, PolY, PolZ);
      //G4cout << "polarization" << dynamicGamma.GetPolarization() << G4endl;
      //   }
  

  dynamicGamma.DumpInfo();

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
      G4cout<< "applicability OK" << G4endl;
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

  //  G4cout<< "table OK" << endl;
  
  /*

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

  G4LowEnergyCompton* gammaLowEProcess =
    (G4LowEnergyCompton*) gammaProcess;
  G4ComptonScattering* gammaStdProcess =
    (G4ComptonScattering*) gammaProcess;
  
  
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
  */
  
  // --------- Test the DoIt 
  
  G4cout << "DoIt in " << material->GetName() << G4endl;

  dynamicGamma.SetKineticEnergy(eEnergy);
  dynamicGamma.SetMomentumDirection(initX,initY,initZ);

  for (G4int iter=0; iter<nIterations; iter++)
    {
      
      step->SetStepLength(1*micrometer);
      


      G4cout  <<  "Iteration = "  <<  iter 
	      << "  -  Step Length = " 
	      << step->GetStepLength()/mm << " mm "
	      << G4endl;
      
      gTrack->SetStep(step); 

      G4StepPoint* preStep  = step->GetPreStepPoint();
      G4StepPoint* postStep = step->GetPostStepPoint();
      G4ThreeVector prePosition = preStep->GetPosition();
      G4ThreeVector postPosition = postStep->GetPosition();

      //G4cout << prePosition << "pre step point "<< G4endl;
      //G4cout << postPosition << "post step point "<< G4endl;

      G4ThreeVector polInitial=dynamicGamma.GetPolarization();

      G4cout << polInitial << " Initial Polarization" << G4endl;

      G4VParticleChange* dummy;
      dummy = gammaProcess->PostStepDoIt(*gTrack, *step);
      G4ParticleChange* particleChange = (G4ParticleChange*) dummy;
      

      // Primary physical quantities 


      //      particleChange->DumpInfo();
      
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
      G4double thetaChange = eChange.theta();


      const G4ThreeVector* momChange =particleChange->GetMomentumDirectionChange();

      G4cout << (momChange->x()) << " " <<  (momChange->y()) << " "  << (momChange->z()) << " "  << G4endl;

      G4cout << eChange << "newdir" << G4endl;
      G4double phiChange = eChange.phi();

      G4double xChange = particleChange->GetPositionChange()->x();
      G4double yChange = particleChange->GetPositionChange()->y();
      G4double zChange = particleChange->GetPositionChange()->z();
      
      //G4cout << "Theta " <<  thetaChange << G4endl;
      //G4cout << "Phi " <<  phiChange << G4endl;

      G4cout << "---- Primary after the step ---- " << G4endl;
 
      G4cout << "Position (x,y,z) = " 
      	     << xChange << "  " 
      	     << yChange << "   " 
      	     << zChange << "   " 
      	     << G4endl;

      G4cout << " Initial Energy " << initEnergy/keV << " keV" << G4endl; 
      G4cout << "---- Energy: " << energyChange/MeV << " MeV,  " 
	     << "(px,py,pz): ("
	     << pxChange/keV << ","
	     << pyChange/keV << "," 
	     << pzChange/keV << ") keV"
	     << G4endl;
      /*      G4cout << "---- Energy loss (dE) = " << dedx/keV << " keV" << G4endl;
	      G4cout << "Stopping power (dE/dx)=" << dedxNow << G4endl;
      */


      G4double electronMass = 511.22*keV; // da inserire la definizione

      G4double Ratio = energyChange/ 
	(initEnergy/(1 + (initEnergy*(1-cos(thetaChange))/electronMass))); 
      // testenergy

      //G4cout << Ratio << "RATIO" << G4endl;
      //G4cout << energyChange/keV << "ENERGY (keV)" << G4endl;


      const G4ThreeVector* polChange=particleChange->GetPolarizationChange();


      //G4cout << pxChange/pChange << "X" << G4endl;
      //G4cout << pyChange/pChange << "Y" << G4endl;
      //G4cout << pzChange/pChange << "Z" << G4endl;

      //G4cout << polChange->x() << "pol X" << G4endl;
      //G4cout << polChange->y() << "pol Y" << G4endl;
      //G4cout << polChange->z() << "pol Z" << G4endl;
      //G4cout << polChange->mag() << "pol mag" << G4endl;


      G4double ScalarProduct = (polChange->x())*(pxChange/pChange)+
	(pyChange/pChange)*(polChange->y())+
	(pzChange/pChange)*(polChange->z());
      
      //G4cout << ScalarProduct << "scalar product" << G4endl;

      hgETot->accumulate(energyChange/keV);
      hgP->accumulate(pChange/keV);
      hgTheta->accumulate(thetaChange);
      hgPhi->accumulate(phiChange);
      hgRapp->accumulate(Ratio);
      
      // Secondaries 

      ntuple1->column("eprimary", initEnergy/keV);
      ntuple1->column("energyf", energyChange/keV);
      ntuple1->column("de", dedx/keV);
      ntuple1->column("dedx", dedxNow/keV);
      ntuple1->column("pxch", pxChange);
      ntuple1->column("pych", pyChange);
      ntuple1->column("pzch", pzChange);
      ntuple1->column("pch", pChange);
      ntuple1->column("polx",(polInitial.x()));
      ntuple1->column("poly",(polInitial.y()));
      ntuple1->column("polz",(polInitial.z()));
      ntuple1->column("polchx",(polChange->x()));
      ntuple1->column("polchy",(polChange->y()));
      ntuple1->column("polchz",(polChange->z()));
      ntuple1->column("thetach", thetaChange);  
      ntuple1->column("phich", phiChange);  
      ntuple1->dumpData(); 

      // Secondaries physical quantities 
      
      hNSec->accumulate(particleChange->GetNumberOfSecondaries());
      hDebug->accumulate(particleChange->GetLocalEnergyDeposit());
      
      G4cout << " secondaries " << 
	particleChange->GetNumberOfSecondaries() << G4endl;

      G4double Etotal = 0.;
      Etotal += energyChange;

      //G4cout << " Total energy" << Etotal << G4endl;      

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
	  
	  if (eKin > initEnergy)
	    {
	      G4cout << "WARNING: eKinFinal > eKinInit " << G4endl;
	      //	     << e
	      //		     << " > " initEnergy 
	      
	    }
	  
	  G4String particleName = 
	    finalParticle->GetDefinition()->GetParticleName();
	  G4cout  << "==== Final " 
		  <<  particleName  <<  " "  
		  << "energy: " <<  e/keV  <<  " keV,  " 
		  << "eKin: " <<  eKin/keV  <<  " keV, " 
		  << "(px,py,pz): ("
		  <<  px/keV  <<  "," 
		  <<  py/keV  <<  ","
		  <<  pz/keV  << ") keV "
		  <<  G4endl;   

	  //	  G4cout << " energia secondaria" << e << G4endl;      
	  
	  heETot->accumulate(eKin/keV);
	  heP->accumulate(p/keV);
          
	  Etotal += eKin;      
	  //G4cout << " energia totale" << Etotal << G4endl;      
	  
	  G4int partType;
	  if (particleName == "e-") partType = 1;
	  else if (particleName == "e+") partType = 2;
	  else if (particleName == "gamma") partType = 3;
	  
	  // Fill the secondaries ntuple

          ntuple2->column("event",iter);
          ntuple2->column("eprimary",initEnergy/keV);
	  ntuple2->column("px", px);
	  ntuple2->column("py", py);
	  ntuple2->column("pz", pz);
	  ntuple2->column("p", p);
	  ntuple2->column("e", e/keV);
	  ntuple2->column("theta", theta);
	  ntuple2->column("ekin", eKin/keV);
	  ntuple2->column("type", partType);
	  
	  ntuple2->dumpData(); 
	  
	  delete particleChange->GetSecondary(i);
	}

      //      G4cout << Etotal/keV << " E total (keV) " << G4endl;
      hSumE->accumulate(Etotal/keV);
      particleChange->Clear();
      
    } 
  
  
  //  G4cout  << "Iteration number: "  <<  iter << G4endl;
  hbookManager->write();
  delete hbookManager;
  
  delete step;

  G4cout << "END OF THE MAIN PROGRAM" << G4endl;
}

















