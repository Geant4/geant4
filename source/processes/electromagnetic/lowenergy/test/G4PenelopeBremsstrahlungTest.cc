//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geSant4/license                                  *
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
// $Id: G4PenelopeBremsstrahlungTest.cc,v 1.5 2004-06-04 06:27:48 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4PenelopeBremsstrahlungTest.cc
//
//      Author:        Francesco Longo
// 
//      Creation date: 04 january 2001
//
//      Modifications: Luciano Pandola  (27 november 2002)
//                     Adapted in order to test G4PenelopeBremsstrahlung
//                     Minor modification in n-tuple filling
//                     Updated analysis to AIDA 3.0 
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4VContinuousDiscreteProcess.hh"
#include "G4VProcess.hh"
#include "G4ProcessManager.hh"
#include "G4RunManager.hh"

#include "G4PenelopeBremsstrahlung.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4eBremsstrahlung.hh"

#include "G4EnergyLossTables.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"
#include "G4ForceCondition.hh"

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
#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"

#include "G4UnitsTable.hh"
#include "AIDA/IManagedObject.h"

#include <memory>
#include "AIDA/IAnalysisFactory.h"
#include "AIDA/ITreeFactory.h"
#include "AIDA/ITree.h"
#include "AIDA/IHistogramFactory.h"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"
#include "AIDA/IHistogram3D.h"
#include "AIDA/ITupleFactory.h"
#include "AIDA/ITuple.h"


int main()
{

  // Setup

  G4int nIterations = 100000;
  G4int materialId = 3;
  G4int test=0;
  G4int tPart=1;
  //G4cout.setf(std::ios::scientific,std::ios::floatfield );

  // -------------------------------------------------------------------

  // ---- HBOOK initialization

  std::auto_ptr< AIDA::IAnalysisFactory > af( AIDA_createAnalysisFactory() );
  std::auto_ptr< AIDA::ITreeFactory > tf (af->createTreeFactory());
  std::auto_ptr< AIDA::ITree > tree (tf->create("pen_br_test.hbook","hbook",false,true));
  G4cout << "Tree store: " << tree->storeName() << G4endl;
  std::auto_ptr< AIDA::ITupleFactory > tpf (af->createTupleFactory(*tree));
  std::auto_ptr< AIDA::IHistogramFactory > hf (af->createHistogramFactory(*tree));
 
  // ---- primary ntuple ------
  AIDA::ITuple* ntuple1 = tpf->create("1","Primary Ntuple","double eprimary,energyf,de,dedx,pxch,pych,pzch,pch,thetach");
  
  // ---- secondary ntuple ------
  AIDA::ITuple* ntuple2 = tpf->create("2","Secondary Ntuple","double eprimary,px_ga,py_ga,pz_ga,p_ga,e_ga,theta_ga");

  // ---- table ntuple ------
  AIDA::ITuple* ntuple3 = tpf->create("3","Mean Free Path Ntuple","double kinen,mfp");

  //--------- Materials definition ---------

  G4Material* Si  = new G4Material("Silicon",   14., 28.055*g/mole, 2.33*g/cm3);
  G4Material* Fe  = new G4Material("Iron",      26., 55.85*g/mole, 7.87*g/cm3);
  G4Material* Cu  = new G4Material("Copper",    29., 63.55*g/mole, 8.96*g/cm3);
  G4Material*  W  = new G4Material("Tungsten", 74., 183.85*g/mole, 19.30*g/cm3);
  G4Material* Pb  = new G4Material("Lead",      82., 207.19*g/mole, 11.35*g/cm3);
  G4Material*  U  = new G4Material("Uranium", 92., 238.03*g/mole, 18.95*g/cm3);
  G4Material* Al  = new G4Material("Aluminum",13.,26.98*g/mole,2.7*g/cm3);
  G4Material* Au  = new G4Material("Gold"    ,79.,196.97*g/mole,19.3*g/cm3);

  G4Element*   H  = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   O  = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
  G4Element*   C  = new G4Element ("Carbon"  , "C", 6. , 12.00*g/mole);
  G4Element*  Cs  = new G4Element ("Cesium"  , "Cs", 55. , 132.905*g/mole);
  G4Element*   I  = new G4Element ("Iodine"  , "I", 53. , 126.9044*g/mole);

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
  G4cout << "Electrons [1] or Positrons [2] ?" << G4endl;
  G4cin >> tPart;
  if ( !(tPart == 1 || tPart == 2)) G4Exception("Wrong input");


  G4cout << "Test AlongStepDoIt [1] or PostStepDoIt [2] ?" << G4endl;
  G4cin >> test;
  if ( !(test == 1 || test == 2)) G4Exception("Wrong input");



  G4cout << "How many interactions? " << G4endl;
  G4cin >> nIterations;

  if (nIterations <= 0) G4Exception("Wrong input");

  G4double initEnergy = 1*MeV; 
  G4double initX = 0.; 
  G4double initY = 0.; 
  G4double initZ = 1.;
  
  G4cout << "Enter the initial particle energy E (MeV)" << G4endl; 
  G4cin >> initEnergy ;


  initEnergy = initEnergy*MeV;
  
  if (initEnergy  <= 0.) G4Exception("Wrong input");

  static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  G4int nMaterials = G4Material::GetNumberOfMaterials();

  G4cout << "Available materials are: " << G4endl;
  for (G4int mat = 0; mat < nMaterials; mat++)
    {
      G4cout << mat << ") "
	     << (*theMaterialTable)[mat]->GetName()
	     << G4endl;
    }
  
  G4cout << "Which material? " << G4endl;
  G4cin >> materialId;
  
  G4Material* material = (*theMaterialTable)[materialId] ;

  G4cout << "The selected material is: "
	 << material->GetName()
	 << G4endl;
  
  G4double dimX = 1*mm;
  G4double dimY = 1*mm;
  G4double dimZ = 1*mm;
  
  // Geometry 
  
  G4Box* theFrame = new G4Box ("Frame",dimX, dimY, dimZ);
  
  G4LogicalVolume* logicalFrame = new G4LogicalVolume(theFrame,
						      (*theMaterialTable)[materialId],
						      "LFrame", 0, 0, 0);
  logicalFrame->SetMaterial(material); 
  
  G4PVPlacement* physicalFrame = new G4PVPlacement(0,G4ThreeVector(),
						   "PFrame",logicalFrame,0,false,0);
  G4RunManager* rm = new G4RunManager();
  G4cout << "World is defined " << G4endl;
  rm->GeometryHasBeenModified();
  rm->DefineWorldVolume(physicalFrame);
  // Particle definitions
  
  G4ParticleDefinition* gamma = G4Gamma::GammaDefinition();
  G4ParticleDefinition* electron = G4Electron::ElectronDefinition();
  G4ParticleDefinition* positron = G4Positron::PositronDefinition();
  
  
  G4ParticleDefinition* realpt = G4Electron::ElectronDefinition();

  if (tPart == 2)
    {
      realpt = G4Positron::PositronDefinition();
    }
  
  G4ProductionCutsTable* cutsTable = G4ProductionCutsTable::GetProductionCutsTable();
  G4ProductionCuts* cuts = cutsTable->GetDefaultProductionCuts();
  G4double cutG=1*micrometer;
  G4double cutE=1*micrometer;
  cuts->SetProductionCut(cutG, 0); //gammas
  cuts->SetProductionCut(cutE, 1); //electrons
  cuts->SetProductionCut(cutE, 2); //positrons
  G4cout << "Cuts are defined " << G4endl;
 
  //G4Gamma::SetEnergyRange(2.5e-4*MeV,1e5*MeV);
  //G4Electron::SetEnergyRange(2.5e-4*MeV,1e5*MeV);
  //G4Positron::SetEnergyRange(2.5e-4*MeV,1e5*MeV);
  
  cutsTable->UpdateCoupleTable();
  //cutsTable->DumpCouples();
  const G4MaterialCutsCouple* theCouple = cutsTable->GetMaterialCutsCouple(material,cuts);

  // Processes 
  
  
  G4int processType;
  G4cout << "Standard [1] or LowEnergy[2] or Penelope [3] Bremsstrahlung?" << G4endl;
  G4cin >> processType;
  if ( !(processType == 1 || processType == 2 || processType == 3))
    {
      G4Exception("Wrong input");
    }

  G4VContinuousDiscreteProcess* bremProcess;

  if (processType == 1)
    {
      bremProcess = new G4eBremsstrahlung();
      G4cout << "The selected model is Standard" << G4endl;
    }
  else if (processType == 2)
    {
      bremProcess = new G4LowEnergyBremsstrahlung();
      G4cout << "The selected model is Low Energy" << G4endl;
    }
  else if (processType == 3)
    {
     
      bremProcess = new G4PenelopeBremsstrahlung();
      G4cout << "The selected model is Penelope" << G4endl;
    }
  
  //----------------
  // process manager  
  //----------------

  // electron or positron
  
  
  G4ProcessManager* ProcessManager = new G4ProcessManager(realpt);
  realpt->SetProcessManager(ProcessManager);
  ProcessManager->AddProcess(bremProcess);
  G4ForceCondition* condition; 

  
  //--------------
  // set ordering   
  //--------------


//   eProcessManager->
//     SetProcessOrdering(theeminusMultipleScattering, idxAlongStep,1);
//   eProcessManager->
//     SetProcessOrdering(theeminusIonisation,         idxAlongStep,2);
      
//   eProcessManager->
//     SetProcessOrdering(theeminusMultipleScattering, idxPostStep,1);
//   eProcessManager->
//     SetProcessOrdering(theeminusIonisation,         idxPostStep,2);
//   eProcessManager->
//     SetProcessOrdering(theeminusBremsstrahlung,     idxPostStep,3);



  // pProcessManager->SetProcessOrderingToFirst(theeplusAnnihilation, idxAtRest);
//   pProcessManager->
//     SetProcessOrdering(theeplusMultipleScattering, idxAlongStep,1);
//   pProcessManager->
//     SetProcessOrdering(theeplusIonisation,         idxAlongStep,2);

//   pProcessManager->
//     SetProcessOrdering(theeplusMultipleScattering, idxPostStep,1);
//   pProcessManager->
//     SetProcessOrdering(theeplusIonisation,         idxPostStep,2);
//   pProcessManager->
//     SetProcessOrdering(theeplusBremsstrahlung,     idxPostStep,3);
//   pProcessManager->
//     SetProcessOrdering(theeplusAnnihilation,       idxPostStep,4);
  
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
  //eEnergy is the KINETIK energy
  G4DynamicParticle dynamicPrimary(realpt,eDirection,eEnergy);

  dynamicPrimary.DumpInfo(0);
  
  // Track 

  G4ThreeVector aPosition(0.,0.,0.);
  G4double aTime = 0. ;
  
  G4Track* eTrack = new G4Track(&dynamicPrimary,aTime,aPosition);
  G4GRSVolume* touche = new G4GRSVolume(physicalFrame, NULL, aPosition);   
  eTrack->SetTouchableHandle(touche); //verificare!!!!!!!!!!!!


  // Step 

  G4Step* step = new G4Step();  
  step->SetTrack(eTrack);

  G4StepPoint* aPoint = new G4StepPoint();
  aPoint->SetPosition(aPosition);
  aPoint->SetMaterial(material);
  aPoint->SetMaterialCutsCouple(theCouple);
  G4double safety = 10000.*cm;
  aPoint->SetSafety(safety);
  step->SetPreStepPoint(aPoint);
  
  // Check applicability
  
  if (! (bremProcess->IsApplicable(*realpt)))
    {
      G4Exception("Not Applicable");
    }
  else 
    {
      G4cout<< "applicability OK" << G4endl;
    }
  
  // Initialize the physics tables (in which material?)
  bremProcess->BuildPhysicsTable(*realpt);

  G4cout<< "table OK" << G4endl;
  
  // Test GetMeanFreePath()
  // E' protected! Il membro accessibile e' DumpMeanFreePath()
  
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
  apttoMaterial = (*theMaterialTable)[materialId] ;
  MaterialName  = apttoMaterial->GetName() ;
  logicalFrame->SetMaterial(apttoMaterial); 
  
  eTrack->SetStep(step);

  G4PenelopeBremsstrahlung* LowEProcess =
    (G4PenelopeBremsstrahlung*) bremProcess;
  G4LowEnergyBremsstrahlung* LowEProcess2 =
    (G4LowEnergyBremsstrahlung*) bremProcess;
  G4eBremsstrahlung* StdProcess =
    (G4eBremsstrahlung*) bremProcess;
 
  
  for (G4int i=0 ; i<pntNum; i++)
    {
      dynamicPrimary.SetKineticEnergy(Tkin[i]);
      if (processType == 3)
	{
	  //meanFreePath=LowEProcess
	  // ->DumpMeanFreePath(*eTrack, sti, condition);
	}
      else if (processType == 2)
	{
	  //meanFreePath=electronLowEProcess2
	  // ->DumpMeanFreePath(*eTrack, sti, condition);
	  
	}
      else if (processType == 1)
	{ 
	  //meanFreePath=StdProcess
	  //  ->GetMeanFreePath(*eTrack, sti, condition); 
	}

      ntuple3->fill(ntuple3->findColumn("kinen"),log10(Tkin[i]));
      ntuple3->fill(ntuple3->findColumn("mfp"),meanFreePath/cm);
      ntuple3->addRow();

    
      //G4cout << Tkin[i]/MeV << " " <<  meanFreePath/cm << G4endl;

    }
  G4cout << "Mean Free Path OK" << G4endl;
  
  // --------- Test the DoIt 
  
  G4cout << "DoIt in " << material->GetName() << G4endl;


  dynamicPrimary.SetKineticEnergy(eEnergy);
  G4int iter;
  for (iter=0; iter<nIterations; iter++)
    {
      
      step->SetStepLength(1*micrometer);
      
      G4cout  <<  "Iteration = "  <<  iter 
	      << "  -  Step Length = " 
	      << step->GetStepLength()/mm << " mm "
	      << G4endl;
      
    
      eTrack->SetStep(step); 
     
 
      //      G4cout  <<  "Iteration = "  <<  iter 
      //	      << "  -  Step Length = " 
      //      << step->GetStepLength()/mm << " mm "
      //      << G4endl;
      
      //G4cout << eTrack->GetStep()->GetStepLength()/mm 
      //     << G4endl;
      
      //G4cout << "Prima" << G4endl;
      G4VParticleChange* dummy;
      if (test==1) dummy = bremProcess->AlongStepDoIt(*eTrack, *step);
      if (test==2) dummy = bremProcess->PostStepDoIt(*eTrack,*step);
      //G4cout << "Dopo" << G4endl;

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

      G4double thetaChange = particleChange->GetMomentumChange()->theta();
      thetaChange = thetaChange/deg; //conversion in degrees
      
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
 
       ntuple1->fill(ntuple1->findColumn("eprimary"),initEnergy/MeV);
       ntuple1->fill(ntuple1->findColumn("energyf"),energyChange/MeV);
       ntuple1->fill(ntuple1->findColumn("de"),dedx/MeV);
       ntuple1->fill(ntuple1->findColumn("dedx"),dedxNow/(MeV/cm));
       ntuple1->fill(ntuple1->findColumn("pxch"),pxChange/MeV);
       ntuple1->fill(ntuple1->findColumn("pych"),pyChange/MeV);
       ntuple1->fill(ntuple1->findColumn("pzch"),pzChange/MeV);
       ntuple1->fill(ntuple1->findColumn("pch"),pChange/MeV);
       ntuple1->fill(ntuple1->findColumn("thetach"),thetaChange);
       ntuple1->addRow();

      // Secondaries physical quantities 
           
      // Secondaries 
      G4cout << " secondaries " << 
	particleChange->GetNumberOfSecondaries() << G4endl;
      G4double px_ga,py_ga,pz_ga,p_ga,e_ga,theta_ga,eKin_ga;
      
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
	  theta = theta/deg; //conversion in degrees
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
	  
	  G4int partType;
	  if (particleName == "e-") {
	    partType = 1;
	  }
	  else if (particleName == "gamma") 
	    {
	      partType = 2;
	      px_ga=px;
	      py_ga=py;
	      pz_ga=pz;
	      p_ga=p;
	      e_ga=e;
	      theta_ga=theta;
	    }
	  else if (particleName == "e+") partType = 3;
	  

	  delete particleChange->GetSecondary(i);
	}
      
      	  // Fill the secondaries ntuple

      // Normalize all to the energy of primary
      // for gammas initEnergy=initP
      ntuple2->fill(ntuple2->findColumn("eprimary"),initEnergy/MeV);
      ntuple2->fill(ntuple2->findColumn("px_ga"),px_ga/MeV);
      ntuple2->fill(ntuple2->findColumn("py_ga"),py_ga/MeV);
      ntuple2->fill(ntuple2->findColumn("pz_ga"),pz_ga/MeV);
      ntuple2->fill(ntuple2->findColumn("p_ga"),p_ga/MeV);
      ntuple2->fill(ntuple2->findColumn("e_ga"),e_ga/MeV);
      ntuple2->fill(ntuple2->findColumn("theta_ga"),theta_ga);
      ntuple2->addRow();
      particleChange->Clear();
      
    } 
  
  
  G4cout  << "Iteration number: "  <<  iter << G4endl;
  
  G4cout << "Committing.............." << G4endl;
  tree->commit();
  G4cout << "Closing the tree........" << G4endl;
  tree->close();
  
  delete step;


  G4cout << "END OF THE MAIN PROGRAM" << G4endl;
  return 0;
}

















