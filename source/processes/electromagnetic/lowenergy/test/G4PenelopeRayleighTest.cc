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
// $Id: G4PenelopeRayleighTest.cc,v 1.1 2002-12-06 16:27:00 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4PenelopeRayleighTest.cc
//
//      Author:        Francesco Longo
// 
//      Creation date: 04 january 2001
//
//      Modifications: Luciano Pandola  (27 november 2002)
//                     Adapted in order to test G4PenelopeRayleigh
//                     Minor modification in n-tuple filling
//                     Updated analysis to AIDA 3.0 
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
#include "G4LowEnergyRayleigh.hh"
#include "G4PenelopeRayleigh.hh"
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


G4int main()
{

  // Setup

  G4int nIterations = 100000;
  G4int materialId = 3;

  //G4cout.setf(G4std::ios::scientific,G4std::ios::floatfield );

  // -------------------------------------------------------------------

  // ---- HBOOK initialization

  G4std::auto_ptr< AIDA::IAnalysisFactory > af( AIDA_createAnalysisFactory() );
  G4std::auto_ptr< AIDA::ITreeFactory > tf (af->createTreeFactory());
  G4std::auto_ptr< AIDA::ITree > tree (tf->create("pen_ray_test.hbook","hbook",false,true));
  G4cout << "Tree store: " << tree->storeName() << G4endl;
  G4std::auto_ptr< AIDA::ITupleFactory > tpf (af->createTupleFactory(*tree));
  G4std::auto_ptr< AIDA::IHistogramFactory > hf (af->createHistogramFactory(*tree));
 
  // ---- primary ntuple ------
  AIDA::ITuple* ntuple1 = tpf->create("1","Primary Ntuple","double eprimary,energyf,de,dedx,pxch,pych,pzch,pch,thetach");
  
  // ---- secondary ntuple ------
  //AIDA::ITuple* ntuple2 = tpf->create("2","Secondary Ntuple","double 
  //eprimary,px_el,py_el,pz_el,p_el,e_el,theta_el,ekin_el,px_po,py_po,pz_po,
  //p_po,e_po,theta_po,ekin_po");

  // ---- table ntuple ------
  AIDA::ITuple* ntuple3 = tpf->create("3","Mean Free Path Ntuple","double kinen,mfp");  

  //--------- Materials definition ---------
  G4Element* SiEl = new G4Element ("SiEl","Si",14.,28.055*g/mole);
  G4Element* FeEl = new G4Element ("FeEl","Fe",26.,55.58*g/mole);
  G4Element* CuEl = new G4Element ("CuEl","Cu",29.,63.55*g/mole);
  G4Element* WEl  = new G4Element ("WEl","W",74.,183.85*g/mole);
  G4Element* PbEl = new G4Element ("PbEl","Pb",82.,207.19*g/mole);
  G4Element* UEl  = new G4Element ("UEl","U",92.,238.03*g/mole);
  G4Element*   H  = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   O  = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
  G4Element*   C  = new G4Element ("Carbon"  , "C", 6. , 12.00*g/mole);
  G4Element*  Cs  = new G4Element ("Cesium"  , "Cs", 55. , 132.905*g/mole);
  G4Element*   I  = new G4Element ("Iodine"  , "I", 53. , 126.9044*g/mole);
  G4Element*   N  = new G4Element ("Nitrogen","N",7.,14.0*g/mole);

  G4Material* Si  = new G4Material("Silicon", 2.33*g/cm3,1);
  Si->AddElement(SiEl,1);
  G4Material* Fe  = new G4Material("Iron",7.87*g/cm3,1);
  Fe->AddElement(FeEl,1);
  G4Material* Cu  = new G4Material("Copper", 8.96*g/cm3,1);
  Cu->AddElement(CuEl,1);
  G4Material*  W  = new G4Material("Tungsten",19.30*g/cm3,1);
  W->AddElement(WEl,1);
  G4Material* Pb  = new G4Material("Lead",11.35*g/cm3,1);
  Pb->AddElement(PbEl,1);
  G4Material*  U  = new G4Material("Uranium",18.95*g/cm3,1);
  U->AddElement(UEl,1);
  G4Material*  maO = new G4Material("Oxygen", 1.1*g/cm3,1);
  maO->AddElement(O,2);
  G4Material* water = new G4Material ("Water" , 1.*g/cm3, 2);
  water->AddElement(H,2);
  water->AddElement(O,1);
  G4Material* ethane = new G4Material ("Ethane" , 0.4241*g/cm3, 2);
  ethane->AddElement(H,6);
  ethane->AddElement(C,2);
  G4Material* csi = new G4Material ("CsI" , 4.53*g/cm3, 2);
  csi->AddElement(Cs,1);
  csi->AddElement(I,1);
  G4Material* Air = new G4Material
    ("Air", 1.2929*kg/m3, 2, kStateGas, 300.00*kelvin, 1.0*atmosphere);
  Air->AddElement(N,0.8);
  Air->AddElement(O,0.2);
 


  // Interactive set-up

  G4cout << "How many interactions? " << G4endl;
  G4cin >> nIterations;

  if (nIterations <= 0) G4Exception("Wrong input");

  G4double initEnergy = 10*keV; 
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
    (*theMaterialTable)[materialId]->GetName()
	<<" is: "<< gamma->GetEnergyCuts()[materialId]/keV << " keV" << G4endl;
  G4cout<<"the cut in energy for e- in: "<<
    (*theMaterialTable)[materialId]->GetName()
	<<" is: "<< electron->GetEnergyCuts()[materialId]/keV << " keV" << G4endl;
  
  // Processes 
  
  
  G4int processType;
  G4cout << "LowEnergy[1] or Penelope [2] Rayleigh?" << G4endl;
  G4cin >> processType;
  if ( !(processType == 1 || processType == 2))
    {
      G4Exception("Wrong input");
    }

  G4VDiscreteProcess* gammaProcess;

  if (processType == 1)
    {
      gammaProcess = new G4LowEnergyRayleigh();
      G4cout << "The selected model is Low Energy" << G4endl;
    }
  else if (processType == 2)
    {
      gammaProcess = new G4PenelopeRayleigh();
      G4cout << "The selected model is Penelope" << G4endl;
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
  G4ForceCondition* condition=0;  //l'ho fissata a zero! E' onesto??

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
  gTrack->SetTouchableHandle(touche); //verificare!!!!!!!!!!!!


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

  G4cout<< "table OK" << G4endl;
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
  
  gTrack->SetStep(step);


  G4LowEnergyRayleigh* gammaLowEProcess2 = (G4LowEnergyRayleigh*) gammaProcess;
  G4PenelopeRayleigh* gammaLowEProcess = (G4PenelopeRayleigh*) gammaProcess;
  
  for (G4int i=0 ; i<pntNum; i++)
    {
      dynamicGamma.SetKineticEnergy(Tkin[i]);
      if (processType == 1)
	{
	  meanFreePath=gammaLowEProcess
	    ->DumpMeanFreePath(*gTrack, sti, condition);
	}
      else if (processType == 2)
	{
	  meanFreePath=gammaLowEProcess2
	    ->DumpMeanFreePath(*gTrack, sti, condition);
	}

      ntuple3->fill(ntuple3->findColumn("kinen"),log10(Tkin[i]));
      ntuple3->fill(ntuple3->findColumn("mfp"),log10(meanFreePath/cm));
      ntuple3->addRow();
      // if (Tkin[i]<10*MeV)
// 	G4cout << "Mean free path @" << Tkin[i] << ": " << meanFreePath/cm << G4endl;
    }
  G4cout << "Mean Free Path OK" << G4endl;
  
  // --------- Test the DoIt 
  
  G4cout << "DoIt in " << material->GetName() << G4endl;


  dynamicGamma.SetKineticEnergy(eEnergy);
  G4int iter;
  for (iter=0; iter<nIterations; iter++)
    {
      
      step->SetStepLength(1*micrometer);
      
      G4cout  <<  "Iteration = "  <<  iter 
	      << "  -  Step Length = " 
	      << step->GetStepLength()/mm << " mm "
	      << G4endl;
      
    
      gTrack->SetStep(step); 
     
      
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
      G4double thetaChange = eChange.theta();
      thetaChange = thetaChange/deg; //conversion in degrees
      G4double xChange = particleChange->GetPositionChange()->x();
      G4double yChange = particleChange->GetPositionChange()->y();
      G4double zChange = particleChange->GetPositionChange()->z();
      //G4double thetaChange = particleChange->GetPositionChange()->theta();
      
      G4cout << "---- Primary after the step ---- " << G4endl;
 
            G4cout << "Position (x,y,z) = " 
      	     << xChange << "  " 
      	     << yChange << "   " 
      	     << zChange << "   " 
      	     << G4endl;

      G4cout << "---- Energy: " << energyChange/MeV << " MeV,  " 
	     << "(px,py,pz): ("
	     << pxChange/MeV << ","
	     << pyChange/MeV << "," 
	     << pzChange/MeV << ") MeV"
	     << G4endl;
      
      G4cout << "---- Energy loss (dE) = " << dedx/keV << " keV" << G4endl;
      //      G4cout << "Stopping power (dE/dx)=" << dedxNow << G4endl;
 
       ntuple1->fill(ntuple1->findColumn("eprimary"),initEnergy/MeV);
       ntuple1->fill(ntuple1->findColumn("energyf"),energyChange/initEnergy);
       ntuple1->fill(ntuple1->findColumn("de"),dedx/MeV);
       ntuple1->fill(ntuple1->findColumn("dedx"),dedxNow/(MeV/cm));
       ntuple1->fill(ntuple1->findColumn("pxch"),pxChange/initEnergy);
       ntuple1->fill(ntuple1->findColumn("pych"),pyChange/initEnergy);
       ntuple1->fill(ntuple1->findColumn("pzch"),pzChange/initEnergy);
       ntuple1->fill(ntuple1->findColumn("pch"),pChange/initEnergy);
       ntuple1->fill(ntuple1->findColumn("thetach"),thetaChange);
       ntuple1->addRow();
     
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

















