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
// $Id: G4ComptonTest.cc,v 1.21 2002-12-13 12:26:00 gcosmo Exp $
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
//      28 Nov 2002  AP  update to AIDA 3 
//      14 Sep 2001  AP  Moved histograms to Lizard 
//      16 Sep 2001  AP  Moved ntuples to Lizard 
//
// -------------------------------------------------------------------
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
#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4Material.hh"
#include "G4VContinuousDiscreteProcess.hh"
#include "G4ProcessManager.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4eBremsstrahlung.hh"
#include "G4LowEnergyCompton.hh"
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

// New Histogramming (from AIDA and Anaphe):
#include <memory> // for the auto_ptr(T>

#include "AIDA/AIDA.h"

int main()
{

  // Setup

  G4cout.setf( ios::scientific, ios::floatfield );

  // -------------------------------------------------------------------

  // ---- HBOOK initialization

  // Creating the analysis factory
  G4std::auto_ptr< AIDA::IAnalysisFactory > af( AIDA_createAnalysisFactory() );

  // Creating the tree factory
  G4std::auto_ptr< AIDA::ITreeFactory > tf( af->createTreeFactory() );

  // Creating a tree mapped to a new hbook file.
  bool readOnly = false;
  bool createFile = true;
  G4std::auto_ptr< AIDA::ITree > tree( tf->create( "comptonhisto.hbook", "hbook", readOnly, createFile ) );
  G4std::cout << "Tree store : " << tree->storeName() << G4std::endl;


  // Next create the nTuples using the factory and open it for writing
  // Creating a tuple factory, whose tuples will be handled by the tree
  G4std::auto_ptr< AIDA::ITupleFactory > tpf( af->createTupleFactory( *tree ) );

  // ---- primary ntuple ------
  // If using Anaphe HBOOK implementation, there is a limitation on the length of the
  // variable names in a ntuple
  AIDA::ITuple* ntuple1 = tpf->create( "1", "Primary tuple", 
			     "float initen, eChng, dedx, dedxNow, pxChng, pyChng, pzChng, pChng, thChng, nElec, nPos, nPhot" );


  // ---- secondary ntuple ------   
  AIDA::ITuple* ntuple2 = tpf->create( "2", "Secondary tuple", 
			     "float px,py,pz,p,e,eKin,theta,phi,partTyp" );


  // ---- secondaries histos ----
  // Creating a histogram factory, whose histograms will be handled by the tree
  G4std::auto_ptr< AIDA::IHistogramFactory > hf( af->createHistogramFactory( *tree ) );

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

  G4Material* Air = new G4Material("Air"  ,  1.290*mg/cm3, 2);
  Air->AddElement(N,0.7);
  Air->AddElement(O,0.3);

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
      G4cout << mat << ") "
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

  G4double initX = 0. * mm; 
  G4double initY = 0. * mm; 
  G4double initZ = 1. * mm;

  G4double dimX = 1 * mm;
  G4double dimY = 1 * mm;
  G4double dimZ = 1 * mm;
  
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
  G4cout 
    << "LowEnergy [1] or Standard [2] or LowEnergyPolarized [3]?" 
    << G4endl;
  cin >> processType;
  if (processType <1 || processType >3 )
    {
      G4Exception("Wrong input");
    }

  G4VContinuousDiscreteProcess* bremProcess;
  G4VContinuousDiscreteProcess* ioniProcess;

  if (processType == 1 || processType == 3)
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
    
  G4ProcessManager* positronProcessManager = new G4ProcessManager(positron);
  positron->SetProcessManager(positronProcessManager);
  positronProcessManager->AddProcess(bremProcess);
  
  // Initialize the physics tables 
  bremProcess->BuildPhysicsTable(*electron);
  ioniProcess->BuildPhysicsTable(*electron);

  // Photon process 
  G4VDiscreteProcess* photonProcess;
  if (processType == 1)
    {
      photonProcess = new G4LowEnergyCompton;
    }
  else if (processType == 2)
    {
      photonProcess = new G4ComptonScattering;
    }
  else
    {
      photonProcess = new G4LowEnergyPolarizedCompton;
    }

  G4ProcessManager* gProcessManager = new G4ProcessManager(gamma);
  gamma->SetProcessManager(gProcessManager);
  gProcessManager->AddProcess(photonProcess);
  photonProcess->BuildPhysicsTable(*gamma);

  // Create a DynamicParticle  
  
  G4double gEnergy = initEnergy*MeV;
  G4ParticleMomentum gDirection(initX,initY,initZ);
  G4DynamicParticle dynamicPhoton(G4Gamma::Gamma(),gDirection,gEnergy);

  // Track 

  G4ThreeVector position(0.,0.,0.);
  G4double time = 0. ;
  G4Track* gTrack = new G4Track(&dynamicPhoton,time,position);

  // Do I really need this?
  G4GRSVolume* touche = new G4GRSVolume(physicalFrame, 0, position);   
  // gTrack->SetTouchable(touche);
 
 // Step 

  G4Step* step = new G4Step();  
  step->SetTrack(gTrack);

  G4StepPoint* point = new G4StepPoint();
  point->SetPosition(position);
  point->SetMaterial(material);
  G4double safety = 10000.*cm;
  point->SetSafety(safety);
  step->SetPreStepPoint(point);

  G4StepPoint* newPoint = new G4StepPoint();
  G4ThreeVector newPosition(0.,0.,1.*mm);
  newPoint->SetPosition(newPosition);
  newPoint->SetMaterial(material);
  step->SetPostStepPoint(newPoint);
  
  // Check applicability
  
  if (! (photonProcess->IsApplicable(*gamma)))
    {
      G4Exception("Not Applicable");
    }

  // --------- Test the DoIt 

  G4cout << "DoIt in material " << material->GetName() << G4endl;

  for (G4int iter=0; iter<nIterations; iter++)
    {
      step->SetStepLength(1*micrometer);

      gTrack->SetStep(step); 
 
      G4cout  <<  "Iteration = "  
	      <<  iter 
	      << "  -  Step Length = " 
	      << step->GetStepLength()/mm 
	      << " mm "
	      << G4endl;

      G4ParticleChange* particleChange = (G4ParticleChange*) photonProcess->PostStepDoIt(*gTrack,*step);
     
      // Primary physical quantities 

      // Primary physical quantities 

      G4double energyChange = particleChange->GetEnergyChange();
      G4double dedx = initEnergy - energyChange ;
      G4double dedxNow = dedx / (step->GetStepLength());
      
      G4ThreeVector eChange = particleChange->CalcMomentum(energyChange,
							   (*particleChange->GetMomentumChange()),
							   particleChange->GetMassChange());

      G4double pxChange = eChange.x();
      G4double pyChange = eChange.y();
      G4double pzChange = eChange.z();
      G4double pChange = sqrt(pxChange*pxChange + pyChange*pyChange + pzChange*pzChange);

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
	  G4double p = sqrt(px*px+py*py+pz*pz);

	  if (eKin > initEnergy)
	    {
	      G4cout << "WARNING: eFinal > eInit " << G4endl;
	    }

	  G4String particleName = finalParticle->GetDefinition()->GetParticleName();
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
	  
	  G4int partType = 0;
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
	    }
	
	  // Fill the secondaries ntuple
          ntuple2->fill( ntuple2->findColumn( "px     " ), px       );
          ntuple2->fill( ntuple2->findColumn( "py     " ), py       );
          ntuple2->fill( ntuple2->findColumn( "pz     " ), pz       );
          ntuple2->fill( ntuple2->findColumn( "p      " ), p        );
          ntuple2->fill( ntuple2->findColumn( "e      " ), e        );
          ntuple2->fill( ntuple2->findColumn( "eKin   " ), eKin     );
          ntuple2->fill( ntuple2->findColumn( "theta  " ), theta    );
          ntuple2->fill( ntuple2->findColumn( "phi    " ), phi      );
          ntuple2->fill( ntuple2->findColumn( "partTyp" ), partType );
	  
          // NEW: Values of attributes are prepared; store them to the nTuple:
          ntuple2->addRow(); // check for returning true ...

	  delete particleChange->GetSecondary(i);
	}

      // Fill the primaries ntuple
      ntuple1->fill( ntuple1->findColumn( "initen"  ), initialEnergy );
      //  ntuple1->fill( ntuple1->findColumn( "eChng"   ), energyChange );
      //  ntuple1->fill( ntuple1->findColumn( "dedx"    ), dedx          );
      //  ntuple1->fill( ntuple1->findColumn( "dedxNow" ), dedxNow       );
      //  ntuple1->fill( ntuple1->findColumn( "pxChng"  ), pxChange      );
      //   ntuple1->fill( ntuple1->findColumn( "pyChng"  ), pyChange      );
      //   ntuple1->fill( ntuple1->findColumn( "pzChng"  ), pzChange      );
      //   ntuple1->fill( ntuple1->findColumn( "pChng"   ), pChange       );
      //   ntuple1->fill( ntuple1->findColumn( "thChng"  ), thetaChange   );
      //   ntuple1->fill( ntuple1->findColumn( "nElec"   ), nElectrons    );
      //   ntuple1->fill( ntuple1->findColumn( "nPos"    ), nPositrons    );
      //   ntuple1->fill( ntuple1->findColumn( "nPhot"   ), nPhotons      );

      //NEW: Values of attributes are prepared; store them to the nTuple:
      ntuple1->addRow();
	          
      particleChange->Clear();
      
    } 
  
  G4cout  << "-----------------------------------------------------"  
	  <<  G4endl;


  // Committing the transaction with the tree
  G4std::cout << "Committing..." << G4std::endl;
  tree->commit();
  G4std::cout << "Closing the tree..." << G4std::endl;
  tree->close();

  delete touche;
  delete step;
  delete gTrack;

  G4cout << "END OF TEST" << G4endl;
}
