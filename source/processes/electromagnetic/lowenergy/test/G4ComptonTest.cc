//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4ComptonTest.cc,v 1.28 2008-04-24 14:14:25 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4ComptonTest
//
//      Author:        Maria Grazia Pia 
//                     Andreas Pfeiffer
// 
//      Creation date: 2 May 2001
//
//      Modifications: 
//      08 Mar 2008  MGP Updated to recent Geant4 interface changes
//      28 Nov 2002  AP  update to AIDA 3 
//      14 Sep 2001  AP  Moved histograms to Lizard 
//      16 Sep 2001  AP  Moved ntuples to Lizard 
//      24 Apr 2008  MGP Upgrate to treat couples correctly from Luciano Pandola's PenelopeComptonTest
//
// -------------------------------------------------------------------
//
// (MGP) The following is obsolete and should be replaced by iAIDA instructions
//
// from: geant4/source/processes/electromagnetic/lowenergy/test/
//
// execute the following lines _before_ gmake, 
// source /afs/cern.ch/sw/lhcxx/share/LHCXX/latest/scripts/setupAnaphe
//
// or, for [t]csh fans:
//
// source /afs/cern.ch/sw/lhcxx/share/LHCXX/latest/scripts/setupAnaphe.csh
//
// both assume that you have the correct PATH to the compiler
//
// [gmake and run your simulation]
//
// to start Lizard:
// /afs/cern.ch/sw/lhcxx/share/LHCXX/latest/scripts/lizard
//
// see also: http://cern.ch/Anaphe 
//
// ********************************************************************

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>

#include "G4Material.hh"
#include "G4VContinuousDiscreteProcess.hh"
#include "G4ProcessManager.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4eBremsstrahlung.hh"
#include "G4LowEnergyCompton.hh"
#include "G4PenelopeCompton.hh"
#include "G4LowEnergyPolarizedCompton.hh"
#include "G4ComptonScattering.hh"
#include "G4LowEnergyIonisation.hh"
#include "G4eIonisation.hh"
#include "G4VeLowEnergyLoss.hh"
#include "G4EnergyLossTables.hh"
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
#include "G4UnitsTable.hh"
#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4RunManager.hh"
#include "G4RegionStore.hh"
#include "G4StateManager.hh"
#include "G4ApplicationState.hh"

// New Histogramming (from AIDA and Anaphe):
#include <memory> // for the auto_ptr(T>

#include "AIDA/AIDA.h"

int main()
{
  // G4String fileName;
  //  G4cout << "Enter histogram file name" << G4endl;
  // G4cin >> fileName;


  // Setup

  //  G4cout.setf( ios::scientific, ios::floatfield );

  // -------------------------------------------------------------------

  // ---- HBOOK initialization

  // Creating the analysis factory
  std::auto_ptr< AIDA::IAnalysisFactory > af( AIDA_createAnalysisFactory() );

  // Creating the tree factory
  std::auto_ptr< AIDA::ITreeFactory > tf( af->createTreeFactory() );

  // Creating a tree mapped to a new hbook file.
  bool readOnly = false;
  bool createFile = true;
  std::auto_ptr< AIDA::ITree > tree( tf->create( "comptonhisto.hbook", "hbook", readOnly, createFile ) );
  std::cout << "Tree store : " << tree->storeName() << std::endl;

  // Next create the nTuples using the factory and open it for writing
  // Creating a tuple factory, whose tuples will be handled by the tree
  std::auto_ptr< AIDA::ITupleFactory > tpf( af->createTupleFactory( *tree ) );

 
  // ---- Primary ntuple ------
  // If using Anaphe HBOOK implementation, there is a limitation on the length of the
  // variable names in a ntuple
  AIDA::ITuple* ntuple1 = tpf->create( "1", "Primary tuple", 
			     "float e0, e1, dedx, dedxNow, px1, py1, pz1, p1, theta1, neminus, neplus, nphoton" );

  // ---- Secondary ntuple ------   
  AIDA::ITuple* ntuple2 = tpf->create( "2", "Secondary tuple", 
			     "float px, py, pz, p, e, ekin, theta, phi, type" );

  // ---- Secondaries histos ----
  // Creating a histogram factory, whose histograms will be handled by the tree
  std::auto_ptr< AIDA::IHistogramFactory > hf( af->createHistogramFactory( *tree ) );

  // Creating an 1-dimensional histogram in the root directory of the tree

  AIDA::IHistogram1D* hEKin;
  hEKin = hf->createHistogram1D("10","Kinetic Energy", 100,0.,10.);
  
  AIDA::IHistogram1D* hP;
  hP = hf->createHistogram1D("20","Momentum", 100,0.,10.);
  
  AIDA::IHistogram1D* hNSec;
  hNSec = hf->createHistogram1D("30","Number of secondaries", 10,0.,10.);
  
  AIDA::IHistogram1D* hDeposit;
  hDeposit = hf->createHistogram1D("40","Local energy deposit", 100,0.,10.);
 
  AIDA::IHistogram1D* hTheta;
  hTheta = hf->createHistogram1D("50","Theta", 100,0.,pi);

  AIDA::IHistogram1D* hPhi;
  hPhi = hf->createHistogram1D("60","Phi", 100,-pi,pi);

  AIDA::IHistogram1D* hE1;
  hE1 = hf->createHistogram1D("70","Scattered particle energy", 100,0.,10.);

  AIDA::IHistogram1D* hEdiff;
  hEdiff = hf->createHistogram1D("80","Energy difference initial-scattered particle", 100,0.,10.);


  // end NEW
  // ================================================================================

  // ==================== end of Histogram and NTuple handling

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
  G4Element*   N  = new G4Element("Nitrogen",   "N" , 7., 14.01*g/mole);

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

  G4Material* air = new G4Material("Air"  ,  1.290*mg/cm3, 2);
  air->AddElement(N,0.7);
  air->AddElement(O,0.3);


  // Interactive set-up

  G4cout << "How many interactions? " << G4endl;
  G4int nIterations;
  G4cin >> nIterations;
  if (nIterations <= 0) G4Exception("Wrong input");

  G4cout << "Enter the initial particle energy E (MeV)" << G4endl; 
  G4double initEnergy; 
  G4cin >> initEnergy ;
  initEnergy = initEnergy * MeV;
  G4double initialEnergy = initEnergy;
  if (initEnergy  <= 0.) G4Exception("Wrong input");



  // Dump the material table
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();  
  G4int nMaterials = G4Material::GetNumberOfMaterials();
  G4cout << "Available materials are: " << G4endl;
  for (G4int mat = 0; mat < nMaterials; mat++)
    {
      G4cout << mat << ". "
	     << (*theMaterialTable)[mat]->GetName()
	     << G4endl;
    }

  G4cout << "Which material? " << G4endl;
  G4int materialId;
  G4cin >> materialId;

  G4Material* material = (*theMaterialTable)[materialId] ;

  G4cout << "The selected material is: "
	 << material->GetName()
	 << G4endl;
  // Geometry 

  G4double dimX = 1 * mm;
  G4double dimY = 1 * mm;
  G4double dimZ = 1 * mm;
  
  G4Box* theFrame = new G4Box ("Frame",dimX, dimY, dimZ);
  G4LogicalVolume* logicalFrame = new G4LogicalVolume(theFrame,material,
						      "LFrame", 0, 0, 0);
  logicalFrame->SetMaterial(material); 
  G4PVPlacement* physicalFrame = new G4PVPlacement(0,G4ThreeVector(),
						   "PFrame",logicalFrame,0,false,0);

  // RunManager 
  G4RunManager* rm = new G4RunManager();

  rm->GeometryHasBeenModified();
  rm->DefineWorldVolume(physicalFrame);
    
  G4cout << "[OK] World is defined " << G4endl;
  G4StateManager::GetStateManager()->SetNewState(G4State_Idle);
  rm->DumpRegion("DefaultRegionForTheWorld"); //this forces the region update!

  // Particle definitions

  G4ParticleDefinition* gamma = G4Gamma::GammaDefinition();
  G4ParticleDefinition* electron = G4Electron::ElectronDefinition();
  G4ParticleDefinition* positron = G4Positron::PositronDefinition();

  G4ProductionCutsTable* cutsTable = G4ProductionCutsTable::GetProductionCutsTable();
  G4ProductionCuts* cuts = cutsTable->GetDefaultProductionCuts();
  if (cuts == 0) G4cout << " G4ProductionCuts* cuts = 0" << G4endl;

  G4double cutG = 1*micrometer;
  G4double cutE = 1*micrometer;
  cuts->SetProductionCut(cutG, gamma);    // photons
  cuts->SetProductionCut(cutE, electron); // electrons
  cuts->SetProductionCut(cutE, positron); // positrons
  cuts->SetProductionCut(cutG, 0); //gammas
  cuts->SetProductionCut(cutE, 1); //electrons
  cuts->SetProductionCut(cutE, 2); //positrons
  //G4double cutAll = 1.*micrometer;
  //cuts->SetProductionCut(cutAll, -1); //all

  G4cout << "Cuts are defined " << G4endl;
  
  // MGP 8/3/2008 - There is something wrong with the Couple, to be investigated 
  G4MaterialCutsCouple* couple = new G4MaterialCutsCouple(material,cuts);
  
  logicalFrame->SetMaterialCutsCouple(couple);

  G4cout << "Recalculation needed: " << couple->IsRecalcNeeded() << G4endl;
  
  cutsTable->UpdateCoupleTable(physicalFrame);
  cutsTable->PhysicsTableUpdated();
  cutsTable->DumpCouples();

  //couple->SetUseFlag(true);
  //cutsTable->UpdateCoupleTable(world);
  //cutsTable->DumpCouples();
  
  //RunManager
  //G4RunManager* rm = new G4RunManager();
  //rm->GeometryHasBeenModified();
  //G4VPhysicalVolume* world(physicalFrame);
  //rm->DefineWorldVolume(world);
  //G4cout << "[OK] World is defined " << G4endl;


  // Processes 

  G4int processType;
  G4cout 
    << "LowEnergy [1] or Penelope [2] or LowEnergyPolarized [3] or Standard [4]?" 
    << G4endl;
  G4cin >> processType;
  if (processType < 1 || processType > 4 ) G4Exception("Wrong input");

  G4VContinuousDiscreteProcess* bremProcess;
  G4VContinuousDiscreteProcess* ioniProcess;

  if (processType < 4)
    {
      bremProcess = new G4LowEnergyBremsstrahlung;
      ioniProcess = new G4LowEnergyIonisation;
    }
  else
    {
      bremProcess = new G4eBremsstrahlung;
      ioniProcess = new G4eIonisation;
    }
 
 
  G4ProcessManager* eProcessManager = new G4ProcessManager(electron);
  electron->SetProcessManager(eProcessManager);
  eProcessManager->AddProcess(bremProcess);
  eProcessManager->AddProcess(ioniProcess);
    
  G4ProcessManager* positronProcessManager = new G4ProcessManager(positron);
  positron->SetProcessManager(positronProcessManager);
  positronProcessManager->AddProcess(new G4eBremsstrahlung);
  positronProcessManager->AddProcess( new G4eIonisation());
 
  // Initialize the physics tables 
  bremProcess->BuildPhysicsTable(*electron);
  ioniProcess->BuildPhysicsTable(*electron);

  // Photon process 
  G4VDiscreteProcess* photonProcess = 0;
  if (processType == 1)
    {
      photonProcess = new G4LowEnergyCompton;
      G4cout << "G4LowEnergyCompton CREATED" << G4endl;
   }
  if (processType == 2)
    {
      photonProcess = new G4PenelopeCompton;
    }
  if (processType == 3)
    {
      photonProcess = new G4LowEnergyPolarizedCompton;
    }  
  if (processType == 4)
    {
      photonProcess = new G4ComptonScattering;
    }

  G4ProcessManager* gProcessManager = new G4ProcessManager(gamma);
  gamma->SetProcessManager(gProcessManager);
  gProcessManager->AddProcess(photonProcess);
  photonProcess->BuildPhysicsTable(*gamma);

  // Primary direction 
  G4double initX = 0. * mm; 
  G4double initY = 0. * mm; 
  G4double initZ = 1. * mm;
  G4ParticleMomentum gDirection(initX,initY,initZ);

  // Check applicability
  
  if (! (photonProcess->IsApplicable(*gamma))) G4Exception("Not Applicable");
 
  // Initialize the physics tables (in which material?)
  photonProcess->BuildPhysicsTable(*gamma);
 

  // --------- Test the DoIt -------------------------------------------------------------------

  G4cout << "DoIt in material " << material->GetName() << G4endl;

  
  G4Track* gTrack = 0;
  G4GRSVolume* touche = 0;
  G4Step* step = 0;
  G4StepPoint* point = 0;
  G4StepPoint* newPoint = 0;
  G4double gEnergy = 0.;
  
  for (G4int iter=0; iter<nIterations; iter++)
    {
      // Primary energy  
      //      gEnergy = initEnergy*MeV*G4UniformRand();
      gEnergy = initEnergy*MeV;
      G4cout << "---- Initial energy = " << gEnergy/MeV << " MeV" << G4endl;

      // Dynamic particle (incident primary)
      G4DynamicParticle dynamicPhoton(G4Gamma::Gamma(),gDirection,gEnergy);

      // Track (incident)
      G4ThreeVector position(0.,0.,0.);
      G4double time = 0. ;
      gTrack = new G4Track(&dynamicPhoton,time,position);

      // Do I really need this?
      touche = new G4GRSVolume(physicalFrame, 0, position);   
      //gTrack->SetTouchable(touche);
  
      // Step 
      step = new G4Step();  
      step->SetTrack(gTrack);
      //const G4MaterialCutsCouple* theCouple(cutsTable->GetMaterialCutsCouple(material, cutsTable->GetDefaultProductionCuts()));
      // if (couple == 0) G4cout << "Couple = 0 in setting Step" <<G4endl;
      //const G4MaterialCutsCouple* couple = new G4MaterialCutsCouple(material,cuts);

      // PreStep point
      point = new G4StepPoint();
      point->SetPosition(position);
      point->SetMaterial(material);
      point->SetMaterialCutsCouple(couple);
      G4double safety = 10000.*cm;
      point->SetSafety(safety);
      step->SetPreStepPoint(point);

      // PostStep point
      newPoint = new G4StepPoint();
      G4ThreeVector newPosition(0.,0.,1.*mm);
      newPoint->SetPosition(newPosition);
      newPoint->SetMaterial(material);
      newPoint->SetMaterialCutsCouple(couple);
      newPoint->SetSafety(safety);
      step->SetPostStepPoint(newPoint);

      // Step length
      step->SetStepLength(1.*micrometer);

      gTrack->SetStep(step); 

      //      const G4MaterialCutsCouple* couple = gTrack->GetMaterialCutsCouple();

      G4cout  <<  "Iteration = "  
	      <<  iter 
	      << "  -  Step Length = " 
	      << step->GetStepLength()/mm 
	      << " mm "
	      << G4endl;

      G4ParticleChange* particleChange = (G4ParticleChange*) photonProcess->PostStepDoIt(*gTrack,*step);
 
      // Primary physical quantities 

      G4double energyChange = particleChange->GetEnergy();
      G4double dedx = gEnergy - energyChange ;
      G4double dedxNow = dedx / (step->GetStepLength());
      
      G4ThreeVector eChange = particleChange->CalcMomentum(energyChange,
							   (*particleChange->GetMomentumDirection()),
							   particleChange->GetMass());
      G4double pxChange = eChange.x();
      G4double pyChange = eChange.y();
      G4double pzChange = eChange.z();
      G4double pChange = std::sqrt(pxChange*pxChange + pyChange*pyChange + pzChange*pzChange);

      G4double xChange = particleChange->GetPosition()->x();
      G4double yChange = particleChange->GetPosition()->y();
      G4double zChange = particleChange->GetPosition()->z();
      G4double thetaChange = particleChange->GetMomentumDirection()->theta();

      G4cout << "---- Energy: " 	    
	     << energyChange/MeV << " MeV,  " 
	     << "(px,py,pz): ("
	     << pxChange/MeV << ","
	     << pyChange/MeV << "," 
	     << pzChange/MeV << ") MeV"
	     << G4endl;

      G4cout << "---- Energy loss (dE) = " << dedx/keV << " keV" << G4endl;
  
      // Primary
      
      hNSec->fill(particleChange->GetNumberOfSecondaries());
      hDeposit->fill(particleChange->GetLocalEnergyDeposit());
      hE1->fill(energyChange);
      hEdiff->fill(dedx);

      G4int nElectrons = 0;
      G4int nPositrons = 0;
      G4int nPhotons = 0;
    


      // Secondaries
  
      for (G4int i = 0; i < (particleChange->GetNumberOfSecondaries()); i++) 
	{  
	  G4Track* finalParticle = particleChange->GetSecondary(i) ;
	  
	  G4double e  = finalParticle->GetTotalEnergy();
	  G4double eKin = finalParticle->GetKineticEnergy();
	  G4double px = (finalParticle->GetMomentum()).x();
	  G4double py = (finalParticle->GetMomentum()).y();
	  G4double pz = (finalParticle->GetMomentum()).z();
	  G4double theta = (finalParticle->GetMomentum()).theta();
	  G4double phi = (finalParticle->GetMomentum()).phi();
	  G4double p = std::sqrt(px*px+py*py+pz*pz);

	  if (eKin > gEnergy)
	    {
	      G4cout << "WARNING: eFinal > eInit " << G4endl;
	    }


	  G4String particleName = finalParticle->GetDefinition()->GetParticleName();
          G4ParticleDefinition* def = finalParticle->GetDefinition();
	  G4cout  << "==== Final " 
		  <<  particleName  << " "  
		  << "energy: " <<  e/MeV  << " MeV,  " 
		  << "eKin: " <<  eKin/MeV  << " MeV, " 
		  << "(px,py,pz): ("
		  <<  px/MeV  << "," 
		  <<  py/MeV  << ","
		  <<  pz/MeV  << ") MeV "
		  <<  G4endl;   
	      
	  hEKin->fill(eKin);
	  hP->fill(p);
	  hTheta->fill(theta);
	  hPhi->fill(phi);
	     	  
	  G4double particleType = -1.;
	  if (def == electron) 
	    {
	      particleType = 1.;
	      nElectrons++;
	    } 
	  if (def == positron) 
	    {
	      particleType = 2.;
	      nPositrons++;
	    }
	    
	  if (def == gamma) 
	    {
	      particleType = 3.;
	      nPhotons++;
	    }
	  
	  // Fill the secondaries ntuple
          ntuple2->fill( ntuple2->findColumn( "px" ), px       );
          ntuple2->fill( ntuple2->findColumn( "py" ), py       );
          ntuple2->fill( ntuple2->findColumn( "pz" ), pz       );
          ntuple2->fill( ntuple2->findColumn( "p" ), p        );
          ntuple2->fill( ntuple2->findColumn( "e" ), e        );
          ntuple2->fill( ntuple2->findColumn( "ekin" ), eKin     );
          ntuple2->fill( ntuple2->findColumn( "theta" ), theta    );
          ntuple2->fill( ntuple2->findColumn( "phi" ), phi      );
          ntuple2->fill( ntuple2->findColumn( "type" ), particleType );
	  
          // NEW: Values of attributes are prepared; store them to the nTuple:
          ntuple2->addRow(); // check for returning true ...
	  
	  delete particleChange->GetSecondary(i);
	} // end loop over secondaries

     // Fill the primaries ntuple
      
      ntuple1->fill( ntuple1->findColumn( "e0"  ), gEnergy );
      ntuple1->fill( ntuple1->findColumn( "e1"   ), energyChange );
      ntuple1->fill( ntuple1->findColumn( "dedx"    ), dedx          );
      ntuple1->fill( ntuple1->findColumn( "dedxNow" ), dedxNow       );
      ntuple1->fill( ntuple1->findColumn( "px1"  ), pxChange      );
      ntuple1->fill( ntuple1->findColumn( "py1"  ), pyChange      );
      ntuple1->fill( ntuple1->findColumn( "pz1"  ), pzChange      );
      ntuple1->fill( ntuple1->findColumn( "p1"   ), pChange       );
      ntuple1->fill( ntuple1->findColumn( "theta1"  ), thetaChange   );
      ntuple1->fill( ntuple1->findColumn( "neminus"   ), (G4double) nElectrons    );
      ntuple1->fill( ntuple1->findColumn( "neplus"    ), (G4double) nPositrons    );
      ntuple1->fill( ntuple1->findColumn( "nphoton"   ), (G4double) nPhotons      );

      //NEW: Values of attributes are prepared; store them to the nTuple:
      ntuple1->addRow();
	          
      particleChange->Clear();

      delete touche;
      delete step;
      delete gTrack;
     
       
    } // end loop over events
  
  G4cout  << "-----------------------------------------------------"  << G4endl;

  // Committing the transaction with the tree
  G4cout << "Committing..." << G4endl;
  tree->commit();
  G4cout << "Closing the tree..." << G4endl;
  tree->close();

  G4cout << "END OF TEST" << G4endl;

	  
}
