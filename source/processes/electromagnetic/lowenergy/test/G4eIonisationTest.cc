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
// $Id: G4eIonisationTest.cc,v 1.13 2001-10-25 22:55:25 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4IonisationTest
//
//      Author:        Maria Grazia Pia
// 
//      Creation date: 20 June 2000
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4Material.hh"
#include "G4VContinuousDiscreteProcess.hh"
#include "G4ProcessManager.hh"
#include "G4LowEnergyIonisationVI.hh"
#include "G4eIonisation.hh"
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

// New Histogramming (from AIDA and Anaphe):
#include "Interfaces/IHistoManager.h"
#include "Interfaces/IHistogram1D.h"
#include "Interfaces/IHistogram2D.h"

// For NtupleTag from Anaphe
#include "NtupleTag/LizardNTupleFactory.h"
using namespace Lizard;

int main()
{

  // Setup

  G4int nIterations = 100000;
  G4int materialId = 3;
  G4int test = 0;

  G4cout.setf( ios::scientific, ios::floatfield );

  // -------------------------------------------------------------------

  // ---- HBOOK initialization

  IHistoManager *hbookManager = createIHistoManager();
  assert (hbookManager != 0);
   hbookManager->selectStore("ioni.hbook");

   // Create a nTuple factory:
  NTupleFactory* factory = createNTupleFactory();

  // ---- primary ntuple ------
 
 // ntuple-name is composition of <fileName>:<dirName>:<ntupleID>
  NTuple* ntuple1 = factory->createC( "ioni1.hbook::1" );
  // Check if successful
 assert (ntuple1 != 0);
  
  // ---- secondary ntuple ------
  
 NTuple* ntuple2 = factory->createC( "ioni2.hbook::2" );
  assert (ntuple2 != 0);
   
  // ---- secondaries histos ----
  IHistogram1D* hEKin;
  hEKin = hbookManager->create1D("10","Kinetic Energy", 100,0.,10.);
 
  IHistogram1D* hP;
  hP = hbookManager->create1D("20","Momentum", 100,0.,10.);
 
  IHistogram1D* hNSec;
  hNSec = hbookManager->create1D("30","Number of secondaries", 10,0.,10.);
  
  IHistogram1D* hDebug;
  hDebug = hbookManager->create1D("40","Debug", 100,0.,200.);

  IHistogram1D* hTheta;
 hTheta = hbookManager->create1D("50","Theta", 100,0.,pi);

 IHistogram1D* hPhi;
  hPhi = hbookManager->create1D("60","Phi", 100,-pi,pi);
  
  //  declare and bind "Quantities" to the Ntuple:

  // First tuple ("Primary"):

  Quantity<float> initialEnergy;
  Quantity<float> energyChange;
  Quantity<float> dedx;
  Quantity<float> dedxNow;
  Quantity<float> pxChange;
  Quantity<float> pyChange;
  Quantity<float> pzChange;
  Quantity<float> pChange;
  Quantity<float> nElectrons;
  Quantity<float> nPositrons;
  Quantity<float> nPhotons;
  //  Add and bind the attributes to the first nTuple
  
  if( !( ntuple1->addAndBind( "eprimary"  , initialEnergy) &&
	 ntuple1->addAndBind( "energyf"   , energyChange ) &&
	 ntuple1->addAndBind( "de"        , dedx         ) &&
	 ntuple1->addAndBind( "dedx"      , dedxNow      ) &&
	 ntuple1->addAndBind( "pxch"      , pxChange     ) &&
	 ntuple1->addAndBind( "pych"      , pyChange     ) &&
	 ntuple1->addAndBind( "pzch"      , pzChange     ) &&
	 ntuple1->addAndBind( "pch"       , pChange      ) &&
	 ntuple1->addAndBind( "eminus"    , nElectrons   ) &&
	 ntuple1->addAndBind( "eplus"     , nPositrons   ) &&
	 ntuple1->addAndBind( "nphotons"  , nPhotons     ) ) ) 
    {
      G4cerr << "Error: unable to add attribute to nTuple1." << G4endl;
      // Must be cleaned up properly before any exit.
      delete ntuple1;
      exit(-1);
    }
  
 // Second nTuple ("Secondary"):
  
  Quantity<float> px;
  Quantity<float> py;
  Quantity<float> pz;
  Quantity<float> p;
  Quantity<float> e;
  Quantity<float> eKin;
  Quantity<float> theta;
  Quantity<float> phi;
  Quantity<float> partType;
  
//  Add and bind the attributes to the second nTuple
  //  if( !( ntuple2->addAndBind( "eprimary",initEnergy ) &&
  if( !( ntuple2->addAndBind( "px"	, px        ) &&
	 ntuple2->addAndBind( "py"	, py        ) &&
	 ntuple2->addAndBind( "pz"	, pz        ) &&
	 ntuple2->addAndBind( "p" 	, p         ) &&
	 ntuple2->addAndBind( "e" 	, e         ) &&
	 ntuple2->addAndBind( "ekin"    , eKin      ) &&
	 ntuple2->addAndBind( "theta"   , theta     ) &&
	 ntuple2->addAndBind( "phi"     , phi       ) &&
	 ntuple2->addAndBind( "type"    , partType  )  ) )
    {
      G4cerr << "Error: unable to add attribute to nTuple2" << G4endl;
      // Must be cleaned up properly before any exit.
      delete ntuple2;
      exit(-1);
    }


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

 // Interactive set-up

  G4cout << "Test AlongStepDoIt [1] or PostStepDoIt [2] ?" << G4endl;
  cin >> test;
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

  initEnergy = initEnergy * MeV;

  if (initEnergy  <= 0.) G4Exception("Wrong input");

  // Dump the material table
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
  
  gamma->SetCuts(1e-3*mm);
  electron->SetCuts(1e-3*mm);
  positron->SetCuts(1e-3*mm);
  
  // Processes 

  G4int processType;
  G4cout << "LowEnergy [1] or Standard [2] Ionisation?" << G4endl;
  cin >> processType;
  if ( !(processType == 1 || processType == 2))
    {
      G4Exception("Wrong input");
    }

  G4VContinuousDiscreteProcess* ionisationProcess;

  if (processType == 1)
    {
      ionisationProcess = new G4LowEnergyIonisationVI;
    }
  else
    {
      ionisationProcess = new G4eIonisation;
    }

  G4ProcessManager* eProcessManager = new G4ProcessManager(electron);
  electron->SetProcessManager(eProcessManager);
  eProcessManager->AddProcess(ionisationProcess);
  
  // Create a DynamicParticle  
  
  G4double eEnergy = initEnergy*MeV;
  G4ParticleMomentum eDirection(initX,initY,initZ);
  G4DynamicParticle dynamicElectron(G4Electron::Electron(),eDirection,eEnergy);  

  // Track 

  G4ThreeVector aPosition(0.,0.,0.);
  G4double aTime = 0. ;

  G4Track* eTrack = new G4Track(&dynamicElectron,aTime,aPosition);

  // MGP Check next statement
  G4Track& aTrack = (*eTrack);

  
  // do I really need this?

  G4GRSVolume* touche = new G4GRSVolume(physicalFrame, NULL, aPosition);   
  eTrack->SetTouchable(touche);
 
 // Step 

  G4Step* step = new G4Step();  
  step->SetTrack(eTrack);

  G4StepPoint* aPoint = new G4StepPoint();
  aPoint->SetPosition(aPosition);
  aPoint->SetMaterial(material);
  G4double safety = 10000.*cm;
  aPoint->SetSafety(safety);
  step->SetPreStepPoint(aPoint);
  step->SetPostStepPoint(aPoint);

  // Check applicability
  
  if (! (ionisationProcess->IsApplicable(*electron))) G4Exception("Not Applicable");

  // Initialize the physics tables 

  ionisationProcess->BuildPhysicsTable(*electron);
	
  // --------- Test the DoIt 

  G4cout << "DoIt in material " << material->GetName() << G4endl;

  for (G4int iter=0; iter<nIterations; iter++)
    {
      step->SetStepLength(1*micrometer);

      G4cout  <<  "Iteration = "  <<  iter << G4endl;
      //	      << "  -  Step Length = " 
      //	      << step->GetStepLength()/mm << " mm "
      //	      << G4endl;

      eTrack->SetStep(step); 
 
      //     G4cout << eTrack.GetStep()->GetStepLength()/mm 
      //	     << G4endl;
 
      G4VParticleChange* dummy= 0;
      if (test == 1) dummy = ionisationProcess->AlongStepDoIt(*eTrack, *step);
      if (test == 2) dummy = ionisationProcess->PostStepDoIt(*eTrack, *step);

      G4ParticleChange* particleChange = (G4ParticleChange*) dummy;
      
      // Primary physical quantities 

      energyChange = particleChange->GetEnergyChange();
      dedx = initEnergy - energyChange ;
      dedxNow = dedx / (step->GetStepLength());
      
      G4ThreeVector eChange = particleChange->CalcMomentum(energyChange,
							   (*particleChange->GetMomentumChange()),
							   particleChange->GetMassChange());

      pxChange  = eChange.x();
      pyChange  = eChange.y();
      pzChange  = eChange.z();
      pChange   = sqrt(pxChange*pxChange + pyChange*pyChange + pzChange*pzChange);

      G4double xChange = particleChange->GetPositionChange()->x();
      G4double yChange = particleChange->GetPositionChange()->y();
      G4double zChange = particleChange->GetPositionChange()->z();

      //  G4cout << "---- Primary after the step ---- " << G4endl;
 
      //      G4cout << "Position (x,y,z) = " 
      //	     << xChange << "  " 
      //	     << yChange << "   " 
      //	     << zChange << "   " 
      //	     << G4endl;

      //  G4cout << "---- Energy: " << energyChange/MeV << " MeV,  " 
      //     << "(px,py,pz): ("
      //     << pxChange/MeV << ","
      //     << pyChange/MeV << "," 
      //     << pzChange/MeV << ") MeV"
      //     << G4endl;

      // G4cout << "---- Energy loss (dE) = " << dedx/keV << " keV" << G4endl;
      
      // Primary

      // Secondaries physical quantities 
      
      hNSec->fill(particleChange->GetNumberOfSecondaries());
      hDebug->fill(particleChange->GetLocalEnergyDeposit());

      nElectrons = 0;
      nPositrons = 0;
      nPhotons = 0;

 // Secondaries

      for (G4int i = 0; i < (particleChange->GetNumberOfSecondaries()); i++) 
	{
	  // The following two items should be filled per event, not
	  // per secondary; filled here just for convenience, to avoid
	  // complicated logic to dump ntuple when there are no secondaries
	  
	  G4Track* finalParticle = particleChange->GetSecondary(i) ;
	  
	  e    = finalParticle->GetTotalEnergy();
	  eKin = finalParticle->GetKineticEnergy();
	  px   = (finalParticle->GetMomentum()).x();
	  py   = (finalParticle->GetMomentum()).y();
	  pz   = (finalParticle->GetMomentum()).z();
	  p   = sqrt(px*px+py*py+pz*pz);
          theta = (finalParticle->GetMomentum()).theta();
	  phi = (finalParticle->GetMomentum()).phi();

	  if (eKin > initEnergy)
	    {	    
	      G4cout << "WARNING: eFinal > eInit " << G4endl;
	    }
	  
	  G4String particleName = finalParticle->GetDefinition()->GetParticleName();
	  
	  G4cout  << "==== Final " 
		  <<  particleName  <<  " "  
		  << "energy: " <<  e/MeV  <<  " MeV,  " 
		  << "eKin: " <<  eKin/MeV  <<  " MeV, " 
		  << "(px,py,pz): ("
		  <<  px/MeV  <<  "," 
		  <<  py/MeV  <<  ","
		  <<  pz/MeV  << ") MeV "
		  <<  G4endl;   
	  
	  hEKin->fill(eKin);
	  hP->fill(p);
	  hTheta->fill(theta);
	  hPhi->fill(phi);

	  partType = 0;
	  if (particleName == "e-") 
	    {
	      partType = 1;
	      nElectrons++;
	    }
	  else if (particleName == "e+") 
	    {
	      partType = 2;
	      nPositrons++;
	    }
	  else if (particleName == "gamma") 
	    {
	      partType = 3;
	      nPhotons++;
	      G4cout << "Fluorescence photon: e = " << e/keV << " keV" << G4endl; 
	    }
	  // NEW: Values of attributes are prepared; store them to the nTuple:
          ntuple2->addRow(); // check for returning true ...
  
	  delete particleChange->GetSecondary(i);
	}
      //NEW: Values of attributes are prepared; store them to the nTuple:
      ntuple1->addRow();
      particleChange->Clear();
      
    } 

 G4cout  << "-----------------------------------------------------"  
	  <<  G4endl;

  //-old  hbookManager->write();

  // Tell the manager which histos to store 
  hbookManager->store("10");
  hbookManager->store("20");
  hbookManager->store("30");
  hbookManager->store("40");
 
// the destructor closes the corresponding file
  delete ntuple1;
  delete ntuple2;

  delete hbookManager;
  delete eTrack;
  delete step;
  delete touche;
  delete aPoint;

  cout << "END OF THE MAIN PROGRAM" << G4endl;
}












