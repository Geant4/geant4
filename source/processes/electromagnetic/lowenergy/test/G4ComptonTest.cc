//
// from: geant4/source/processes/electromagnetic/lowenergy/test/
//
// execute the following lines _before_ gmake, 
// select the Anaphe version you want to use (3.6.4-sec for RH61 "old" compiler,
// 3.6.4 for RH61, new compiler (gcc-2.95.2)):
//
// export PATH=$PATH:/afs/cern.ch/sw/lhcxx/specific/redhat61/egcs_1.1.2/3.6.4-sec/bin
// source /afs/cern.ch/sw/lhcxx/share/LHCXX/3.6.4-sec/install/sharedstart.sh
//
// (for the new compiler, use (sh derivates):
// export PATH=$PATH:/afs/cern.ch/sw/lhcxx/specific/redhat61/gcc-2.95.2/3.6.4/bin
// source /afs/cern.ch/sw/lhcxx/share/LHCXX/3.6.4/install/sharedstart.sh
// )
//
// or, for [t]csh fans:
//
// setenv PATH ${PATH}:/afs/cern.ch/sw/lhcxx/specific/redhat61/egcs_1.1.2/3.6.4-sec/bin
// source /afs/cern.ch/sw/lhcxx/share/LHCXX/3.6.4-sec/install/sharedstart.csh
//
// [gmake and run your simulation]
//
// to start Lizard:
// startLizard.sh --noLicense  
//
// see also: http://cern.ch/Anaphe 
//
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
// $Id: G4ComptonTest.cc,v 1.15 2001-11-30 01:26:19 pia Exp $
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
//      14 Sep 2001  AP  Moved histograms to Lizard 
//      16 Sep 2001  AP  Moved ntuples to Lizard 
//
// -------------------------------------------------------------------

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
#include "Interfaces/IHistoManager.h"
#include "Interfaces/IHistogram1D.h"
#include "Interfaces/IHistogram2D.h"

// For NtupleTag from Anaphe
#include "NtupleTag/LizardNTupleFactory.h"
using namespace Lizard;

// Old Histogramming and Ntuple:
//#include "CLHEP/Hist/TupleManager.h"
//#include "CLHEP/Hist/HBookFile.h"
//#include "CLHEP/Hist/Histogram.h"
//#include "CLHEP/Hist/Tuple.h"
//HepTupleManager* hbookManager;

int main()
{

  // Setup

  G4cout.setf( ios::scientific, ios::floatfield );

  // -------------------------------------------------------------------

  // ---- HBOOK initialization

  //-old  hbookManager = new HBookFile("comptontest.hbook", 58);

  IHistoManager *hbookManager = createIHistoManager();
  assert (hbookManager != 0);
  hbookManager->selectStore("comptonhisto.hbook");

  // Create a nTuple factory:
  NTupleFactory* factory = createNTupleFactory();

  // Next create the nTuples using the factory and open it for writing

  // ---- primary ntuple ------
  //-old Tuple* ntuple1 = hbookManager->ntuple("Primary Ntuple");


  // ntuple-name is composition of <fileName>:<dirName>:<ntupleID>
  NTuple* ntuple1 = factory->createC( "comptonhisto.hbook::1" );
  // Check if successful
  assert ( ntuple1 != 0 );

  // ---- secondary ntuple ------   
  //-old HepTuple* ntuple2 = hbookManager->ntuple("Secondary Ntuple");
  NTuple* ntuple2 = factory->createC( "comptonhisto.hbook::2" );
  // Check if successful
  assert ( ntuple2 != 0 );

  // We do not use the factory any longer; it can be deleted.
  delete factory;
  factory = 0;

  // ---- secondaries histos ----
  IHistogram1D* hEKin;
  hEKin = hbookManager->create1D("10","Kinetic Energy", 100,0.,10.);
  
  IHistogram1D* hP;
  hP = hbookManager->create1D("20","Momentum", 100,0.,10.);
  
  IHistogram1D* hNSec;
  hNSec = hbookManager->create1D("30","Number of secondaries", 10,0.,10.);
  
  IHistogram1D* hDeposit;
  hDeposit = hbookManager->create1D("40","Local energy deposit", 100,0.,10.);
 
  IHistogram1D* hTheta;
  hTheta = hbookManager->create1D("50","Theta", 100,0.,pi);

  IHistogram1D* hPhi;
  hPhi = hbookManager->create1D("60","Phi", 100,-pi,pi);

  // ================================================================================
  // NEW: declare and bind "Quantities" to the Ntuple:

  // Declare the quantities which reflect the attributes of the nTuple(s).
  // These quantities are used as variables of their base-types (float),
  // the binding to their nTuple makes sure the updated values are stored.
  // Don't forget to remove the declarations of these variables in the
  // code below.

  // First tuple ("Primary"):
  Quantity<float> initialEnergy;
  Quantity<float> energyChange;
  Quantity<float> dedx;
  Quantity<float> dedxNow;
  Quantity<float> pxChange;
  Quantity<float> pyChange;
  Quantity<float> pzChange;
  Quantity<float> pChange;
  Quantity<float> thetaChange;
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
	 ntuple1->addAndBind( "thetach"   , thetaChange  ) &&
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
  initialEnergy = initEnergy;
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
  gTrack->SetTouchable(touche);
 
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

      //-ap already defined as attributes (see above)
      /*-ap G4double */ energyChange = particleChange->GetEnergyChange();
      /*-ap G4double */ dedx = initEnergy - energyChange ;
      /*-ap G4double */ dedxNow = dedx / (step->GetStepLength());
      
      G4ThreeVector eChange = particleChange->CalcMomentum(energyChange,
							   (*particleChange->GetMomentumChange()),
							   particleChange->GetMassChange());

      /*-ap G4double */ pxChange = eChange.x();
      /*-ap G4double */ pyChange = eChange.y();
      /*-ap G4double */ pzChange = eChange.z();
      /*-ap G4double */ pChange = sqrt(pxChange*pxChange + pyChange*pyChange + pzChange*pzChange);

      G4double xChange = particleChange->GetPositionChange()->x();
      G4double yChange = particleChange->GetPositionChange()->y();
      G4double zChange = particleChange->GetPositionChange()->z();
      /*-ap G4double */ thetaChange = particleChange->GetPositionChange()->theta();

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

//-ap	   ntuple1->column("eprimary", initEnergy);
//-ap	   ntuple1->column("energyf", energyChange);
//-ap	   ntuple1->column("de", dedx);
//-ap	   ntuple1->column("dedx", dedxNow);
//-ap	   ntuple1->column("pxch", pxChange);
//-ap	   ntuple1->column("pych", pyChange);
//-ap	   ntuple1->column("pzch", pzChange);
//-ap	   ntuple1->column("pch", pChange);  
//-ap	   ntuple1->column("thetach", thetaChange);  
      
      hNSec->fill(particleChange->GetNumberOfSecondaries());
      hDeposit->fill(particleChange->GetLocalEnergyDeposit());
      
      /*-ap G4int */ nElectrons = 0;
      /*-ap G4int */ nPositrons = 0;
      /*-ap G4int */ nPhotons = 0;

      // Secondaries
      
      for (G4int i = 0; i < (particleChange->GetNumberOfSecondaries()); i++) 
	{  
	  G4Track* finalParticle = particleChange->GetSecondary(i) ;
	  
	  /*-ap G4double */ e  = finalParticle->GetTotalEnergy();
	  /*-ap G4double */ eKin = finalParticle->GetKineticEnergy();
	  /*-ap G4double */ px = (finalParticle->GetMomentum()).x();
	  /*-ap G4double */ py = (finalParticle->GetMomentum()).y();
	  /*-ap G4double */ pz = (finalParticle->GetMomentum()).z();
	  /*-ap G4double */ theta = (finalParticle->GetMomentum()).theta();
	  /*-ap G4double */ phi = (finalParticle->GetMomentum()).phi();
	  /*-ap G4double */ p = sqrt(px*px+py*py+pz*pz);

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
	  
	  /*-ap G4int */ partType = 0;
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
//-ap	       ntuple2->column("eprimary",initEnergy);
//-ap	       ntuple2->column("px", px);
//-ap	       ntuple2->column("py", py);
//-ap	       ntuple2->column("pz", pz);
//-ap	       ntuple2->column("p", p);
//-ap	       ntuple2->column("e", e);
//-ap	       ntuple2->column("ekin", eKin);
//-ap	       ntuple2->column("theta", theta);
//-ap	       ntuple2->column("phi", phi);
//-ap	       ntuple2->column("type", partType);
//-ap	       
//-ap	       ntuple2->dumpData(); 
	  
          // NEW: Values of attributes are prepared; store them to the nTuple:
          ntuple2->addRow(); // check for returning true ...

	  delete particleChange->GetSecondary(i);
	}

//-ap	   ntuple1->column("nelectrons",nElectrons);
//-ap	   ntuple1->column("npositrons",nPositrons);
//-ap	   ntuple1->column("nphotons",nPhotons);
//-ap	   ntuple1->dumpData(); 

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
  hbookManager->store("50");
  hbookManager->store("60");

  // the destructor closes the corresponding file
  delete ntuple1;
  delete ntuple2;

  // the destructor closes the file
  delete hbookManager;

  delete touche;
  delete step;
  delete gTrack;

  G4cout << "END OF TEST" << G4endl;
}
