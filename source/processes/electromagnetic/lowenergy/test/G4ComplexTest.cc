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
// $Id: G4ComplexTest.cc,v 1.14 2002-07-19 17:34:39 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4ComplexTest
//                     This test provide AlongStepDoIt and PostStepDoIt 
//                     tests for electromagnetic processes. The input
//                     data have to be describe in ASCII file
//
//      Author:        V.Ivanchenko on base of Maria Grazia Pia tests
// 
//      Creation date: 8 May 2001
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

#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyRayleigh.hh"

#include "G4hLowEnergyIonisation.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

#include "G4MultipleScattering.hh"

#include "G4EnergyLossTables.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"
#include "G4AntiProton.hh"
#include "G4Proton.hh"
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

#include "AIDA/IAnalysisFactory.h"

#include "AIDA/ITreeFactory.h"
#include "AIDA/ITree.h"

#include "AIDA/IHistogramFactory.h"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"
//#include "AIDA/IHistogram3D.h"

#include "AIDA/ITupleFactory.h"
#include "AIDA/ITuple.h"
#include "G4Timer.hh"

int main(int argc,char** argv)
{

  // -------------------------------------------------------------------
  // Setup

  G4int  nEvt        = 100;
  G4int  nPart       =-1;
  G4String  nameMat  = "Si";
  G4int  nProcess    = 0;
  G4bool usepaw      = false;
  G4bool postDo      = true;
  G4bool lowE        = true;
  G4int verbose      = 0;
  G4double gEnergy   = 0.1*MeV;
  G4String hFile     = "";
  G4double theStep   = 1.0*micrometer;
  G4double range     = 1.0*micrometer;
  G4double cutG      = 1.0*micrometer;
  G4double cutE      = 1.0*micrometer;
  G4Material* material = 0; 
  G4String name[6] = {"Ionisation", "Bremsstrahlung", "Compton", 
                      "GammaConversion", "PhotoElectric", "Raylaigh"};


  G4cout.setf( ios::scientific, ios::floatfield );

  // -------------------------------------------------------------------
  // Control on input

  if(argc < 2) {
    G4cout << "Input file is not specified! Exit" << G4endl;
    exit(1);
  }

  ifstream* fin = new ifstream();
  string fname = argv[1];
  fin->open(fname.c_str());
  if( !fin->is_open()) {
    G4cout << "Input file <" << fname << "> does not exist! Exit" << G4endl;
    exit(1);
  }

  // -------------------------------------------------------------------
  //--------- Materials definition ---------

  G4Material* m;
  m = new G4Material("Be",    4.,  9.01*g/mole, 1.848*g/cm3);
  m = new G4Material("Graphite",6., 12.00*g/mole, 2.265*g/cm3 );
  m = new G4Material("Al", 13., 26.98*g/mole, 2.7 *g/cm3);
  m = new G4Material("Si",   14., 28.055*g/mole, 2.33*g/cm3);
  m = new G4Material("LAr",   18., 39.95*g/mole, 1.393*g/cm3);
  m = new G4Material("Fe",      26., 55.85*g/mole, 7.87*g/cm3);
  m = new G4Material("Cu",    29., 63.55*g/mole, 8.96*g/cm3);
  m = new G4Material("W", 74., 183.85*g/mole, 19.30*g/cm3);
  m = new G4Material("Pb",      82., 207.19*g/mole, 11.35*g/cm3);
  m = new G4Material("U", 92., 238.03*g/mole, 18.95*g/cm3);

  G4Element*   H  = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   O  = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
  G4Element*   C  = new G4Element ("Carbon"  , "C", 6. , 12.00*g/mole);
  G4Element*  Cs  = new G4Element ("Cesium"  , "Cs", 55. , 132.905*g/mole);
  G4Element*   I  = new G4Element ("Iodide"  , "I", 53. , 126.9044*g/mole);

  m = new G4Material("O2", 8., 16.00*g/mole, 1.1*g/cm3);

  m = new G4Material ("Water" , 1.*g/cm3, 2);
  m->AddElement(H,2);
  m->AddElement(O,1);

  m = new G4Material ("Ethane" , 0.4241*g/cm3, 2);
  m->AddElement(H,6);
  m->AddElement(C,2);
  
  m = new G4Material ("CsI" , 4.53*g/cm3, 2);
  m->AddElement(Cs,1);
  m->AddElement(I,1);

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  G4int nMaterials = G4Material::GetNumberOfMaterials();

  G4cout << "Available materials are: " << G4endl;
  G4int mat;
  for (mat = 0; mat < nMaterials; mat++) {
    G4cout << mat << ") " << (*theMaterialTable)[mat]->GetName() << G4endl;
  }

  G4cout << "Available processes are: " << G4endl;
  for (mat = 0; mat < 6; mat++) {
    G4cout << mat << ") " << name[mat] << G4endl;
  }

  // Particle definitions

  G4ParticleDefinition* gamma = G4Gamma::GammaDefinition();
  G4ParticleDefinition* electron = G4Electron::ElectronDefinition();
  G4ParticleDefinition* positron = G4Positron::PositronDefinition();
  G4ParticleDefinition* proton   = G4Proton::ProtonDefinition();
  G4ParticleDefinition* antiproton = G4AntiProton::AntiProtonDefinition();
  G4ParticleDefinition* part = gamma;

  // Geometry 

  G4double initX = 0.; 
  G4double initY = 0.; 
  G4double initZ = 1.;
  G4double dimX = 100.0*cm;
  G4double dimY = 100.0*cm;
  G4double dimZ = 100.0*cm;  

  G4Box* sFrame = new G4Box ("Box",dimX, dimY, dimZ);
  G4LogicalVolume* lFrame = new G4LogicalVolume(sFrame,material,"Box",0,0,0);
  G4PVPlacement* pFrame = new G4PVPlacement(0,G4ThreeVector(),"Box",
                                            lFrame,0,false,0);

  // -------------------------------------------------------------------
  // ---- Read input file
  G4cout << "Available commands are: " << G4endl;
  G4cout << "#events" << G4endl;
  G4cout << "#particle" << G4endl;
  G4cout << "#energy(MeV)" << G4endl;
  G4cout << "#cutG(mm)" << G4endl;
  G4cout << "#cutE(mm)" << G4endl;
  G4cout << "#range(mm)" << G4endl;
  G4cout << "#step(mm)" << G4endl;
  G4cout << "#material" << G4endl;
  G4cout << "#process" << G4endl;
  G4cout << "#domain" << G4endl;
  G4cout << "#test" << G4endl;
  G4cout << "#paw" << G4endl;
  G4cout << "#verbose" << G4endl;
  G4cout << "#run" << G4endl;
  G4cout << "#exit" << G4endl;
  G4cout << pFrame << G4endl;

  G4ProcessManager *elecManager, *positManager, *gammaManager, 
                   *protManager, *aprotManager;

  elecManager = new G4ProcessManager(electron);
  electron->SetProcessManager(elecManager);

  positManager = new G4ProcessManager(positron);
  positron->SetProcessManager(positManager);

  gammaManager = new G4ProcessManager(gamma);
  gamma->SetProcessManager(gammaManager);

  protManager = new G4ProcessManager(proton);
  proton->SetProcessManager(protManager);

  aprotManager = new G4ProcessManager(antiproton);
  antiproton->SetProcessManager(aprotManager);

  G4eIonisation* elecSTion = 0;
  G4eBremsstrahlung* elecSTbr = 0;
  G4LowEnergyIonisation* elecLEion = 0;
  G4LowEnergyBremsstrahlung* elecLEbr = 0;
  G4bool ionis = true;

  string line, line1;
  G4bool end = true;
  for(G4int run=0; run<100; run++) {
    do {
      (*fin) >> line;
      G4cout << "Next line " << line << G4endl;
      if(line == "#events") {
        (*fin) >> nEvt;
        if(nEvt < 1) nEvt = 1;
      } else if(line == "#particle") {
        (*fin) >> nPart;
      } else if(line == "#energy(MeV)") {
        (*fin) >> gEnergy;
        gEnergy *= MeV;
      } else if(line == "#cutG(mm)") {
        (*fin) >> cutG;
        cutG *= mm;
      } else if(line == "#cutE(mm)") {
        (*fin) >> cutE;
        cutE *= mm;
      } else if(line == "#range(mm)") {
        (*fin) >> range;
        range *= mm;
      } else if(line == "#step(mm)") {
        (*fin) >> theStep;
        theStep *= mm;
      } else if(line == "#material") {
        (*fin) >> nameMat;
      } else if(line == "#process") {
        (*fin) >> nProcess;
      } else if(line == "#domain") {
        (*fin) >> line1;
        if(line1 == "lowenergy") {lowE = true;}
        else                     {lowE = false;}
      } else if(line == "#test") {
        (*fin) >> line1;
        if(line1 == "PostStep") {postDo = true;}
	else                    {postDo = false;}
      } else if(line == "#paw") {
        usepaw = true;
        (*fin) >> hFile;
      } else if(line == "#run") {
        break;
      } else if(line == "#verbose") {
        (*fin) >> verbose;
      } else if(line == "#exit") {
        end = false;
        break;
      }
    } while(end);

    if(!end) break;

    G4cout << "###### Start new run # " << run << "     #####" << G4endl;
    if(nPart == 0) {
      part = gamma;
    } else if(nPart == 1) {
      part = electron;
    } else if(nPart == 2) {
      part = positron;
    } else if(nPart == 3) {
      part = proton;
    } else if(nPart == 4) {
      part = antiproton;
    } else {
      G4cout << "Particle #" << nPart 
             << " is absent in the list of particles: Exit" << G4endl;
      break;
    }
    if(nProcess < 0 || nProcess > 5) {
      G4cout << "Process #" << nProcess 
             << " is absent in the list of processes: Exit" << G4endl;
      break;
    }

    for (mat = 0; mat < nMaterials; mat++) {
      material = (*theMaterialTable)[mat];
      if(nameMat == material->GetName()) break;
    }

    G4cout << "The particle: " << part->GetParticleName() << G4endl;
    G4cout << "The energy:   " << gEnergy/MeV << " MeV" << G4endl;
    G4cout << "The material: " << material->GetName() << G4endl;
    G4cout << "The cut on e-:" << cutE/mm << " mm" << G4endl;
    G4cout << "The cut on g: " << cutG/mm << " mm" << G4endl;
    G4cout << "The step:     " << theStep/mm << " mm" << G4endl;
    if(postDo && lowE) {
      G4cout << "Test of PostStepDoIt  for " << name[nProcess] 
             << " for lowenergy" << G4endl;
    } else if(postDo && !lowE) {
      G4cout << "Test of PostStepDoIt  for " << name[nProcess] 
             << " for standard" << G4endl;
    } else if(!postDo && !lowE) {
      G4cout << "Test of AlongStepDoIt  for " << name[nProcess] 
             << " for standard" << G4endl;
    } else if(!postDo && lowE) {
      G4cout << "Test of AlongStepDoIt  for " << name[nProcess] 
             << " for lowenergy" << G4endl;
    }

    // -------------------------------------------------------------------
    // ---- HBOOK initialization

    // Creating the analysis factory
    G4std::auto_ptr< IAnalysisFactory > af( AIDA_createAnalysisFactory() );

    // Creating the tree factory
    G4std::auto_ptr< ITreeFactory > tf( af->createTreeFactory() );

    // Creating a tree mapped to a new hbook file.
    G4std::auto_ptr< ITree > tree( tf->create( hFile, false, true, "hbook" ) );
    G4cout << "Tree store : " << tree->storeName() << G4endl;
 
    // Creating a tuple factory, whose tuples will be handled by the tree
    //    G4std::auto_ptr< ITupleFactory > tpf( af->createTupleFactory( *tree ) );

    IHistogram1D* hist[4];
    //ITuple* ntuple1 = 0;
    //ITuple* ntuple2 = 0;

    if(usepaw) {

      // ---- primary ntuple ------
      // If using Anaphe HBOOK implementation, there is a limitation on the length of the
      // variable names in a ntuple
      //ntuple1 = tpf->create( "100", "Primary", "float ekin, dedx" );
      //ntuple2 = tpf->create( "101", "Secondary", "float ekin, dedx" );


      // Creating a histogram factory, whose histograms will be handled by the tree
      G4std::auto_ptr< IHistogramFactory > hf( af->createHistogramFactory( *tree ) );

      // Creating an 1-dimensional histogram in the root directory of the tree

      hist[0] = hf->create1D("11","Kinetic Energy (T/T0)", 50,0.,1.0); 
      hist[1] = hf->create1D("12","Momentum (MeV/c)", 50,0.,gEnergy*0.1/MeV);
      hist[2] = hf->create1D("13","Number of secondaries", 20,-0.5,19.5);
      hist[3] = hf->create1D("14","Energy deposition (MeV)", 50,0.,gEnergy*0.1/MeV);

      G4cout<< "Histograms is initialised" << G4endl;
    }

    G4Timer* timer = new G4Timer();
    timer->Start();

    gamma->SetCuts(cutG);
    electron->SetCuts(cutE);
    //    positron->SetCuts(cutE);
  
    // Processes 

    G4VDiscreteProcess*            dProcess;
    G4VContinuousDiscreteProcess* cdProcess;
    dProcess = 0;
    cdProcess = 0;

    G4cout  <<  "Start BuildPhysicsTable"  <<  G4endl;; 

    if(lowE) {
      if(ionis) {
        elecLEion = new G4LowEnergyIonisation();
        elecLEbr = new G4LowEnergyBremsstrahlung();
        elecManager->AddProcess(elecLEion);
        elecManager->AddProcess(elecLEbr);
        elecLEion->BuildPhysicsTable(*electron);
        elecLEbr->BuildPhysicsTable(*electron);
        ionis = false;
      }
      if(nPart == 0) {
        if(nProcess == 2) dProcess =  new G4LowEnergyCompton();
        if(nProcess == 3) dProcess =  new G4LowEnergyGammaConversion();
        if(nProcess == 4) dProcess =  new G4LowEnergyPhotoElectric();
        if(nProcess == 5) dProcess =  new G4LowEnergyRayleigh(); 
        if(dProcess) {
          gammaManager->AddProcess(dProcess);
          dProcess->BuildPhysicsTable(*gamma);
	}
      } else if(nPart == 1) {
        if(nProcess == 0) cdProcess =  elecLEion;
        if(nProcess == 1) cdProcess =  elecLEbr;
      } else if(nPart == 3) {
        if(nProcess == 0) {
          cdProcess = new G4hLowEnergyIonisation();
          protManager->AddProcess(cdProcess);
          cdProcess->BuildPhysicsTable(*proton);
	}
      } else if(nPart == 4) {
        if(nProcess == 0) {
          cdProcess = new G4hLowEnergyIonisation();
          aprotManager->AddProcess(cdProcess);
          cdProcess->BuildPhysicsTable(*antiproton);
	}
      }

    } else {
      if(ionis) {
        elecSTion = new G4eIonisation();
        elecSTbr  = new G4eBremsstrahlung();
        elecManager->AddProcess(elecSTion);
        elecManager->AddProcess(elecSTbr);
        elecSTion->BuildPhysicsTable(*electron);
        elecSTbr->BuildPhysicsTable(*electron);
        ionis = false;
      }
      if(nPart == 0) {
        if(nProcess == 2) dProcess =  new G4ComptonScattering();
        if(nProcess == 3) dProcess =  new G4GammaConversion();
        if(nProcess == 4) dProcess =  new G4PhotoElectricEffect();
        if(dProcess) {
          gammaManager->AddProcess(dProcess);
          dProcess->BuildPhysicsTable(*gamma);
	}
      } else if(nPart == 1) {
        if(nProcess == 0) cdProcess =  elecSTion;
        if(nProcess == 1) cdProcess =  elecSTbr;
      } else if(nPart == 2) {
        G4eIonisation* pSTion = new G4eIonisation();
        G4eBremsstrahlung* pSTbr = new G4eBremsstrahlung();
        if(nProcess == 0) cdProcess =  pSTion;
        if(nProcess == 1) cdProcess =  pSTbr;
      } else if(nPart == 3) {
        if(nProcess == 0) {
          cdProcess = new G4hIonisation();
          protManager->AddProcess(cdProcess);
          cdProcess->BuildPhysicsTable(*proton);
	}
      } else if(nPart == 4) {
        if(nProcess == 0) {
          cdProcess = new G4hIonisation();
          aprotManager->AddProcess(cdProcess);
          cdProcess->BuildPhysicsTable(*antiproton);
	}
      }
    }

    G4cout  <<  "Physics tables are built"  <<  G4endl;; 

    timer->Stop();
    G4cout << "  "  << *timer << G4endl;
    delete timer;

  
    // Control on processes
    if(postDo && !dProcess && !cdProcess) {
        G4cout << "Discret Process is not found out! Exit" << G4endl;
        break;
    }
    if(!postDo && !cdProcess) {
        G4cout << "Continues Discret Process is not found out! Exit" << G4endl;
        break;
    }

    // Create a DynamicParticle  
  
    G4ParticleMomentum gDir(initX,initY,initZ);
    G4DynamicParticle dParticle(part,gDir,gEnergy);

    // Track 
    G4ThreeVector aPosition(0.,0.,0.);
    G4double aTime = 0. ;

    G4Track* gTrack;
    gTrack = new G4Track(&dParticle,aTime,aPosition);

    // Step 

    G4Step* step;
    step = new G4Step();  
    step->SetTrack(gTrack);

    G4StepPoint *aPoint, *bPoint;
    aPoint = new G4StepPoint();
    aPoint->SetPosition(aPosition);
    aPoint->SetMaterial(material);
    G4double safety = 10000.*cm;
    aPoint->SetSafety(safety);
    step->SetPreStepPoint(aPoint);

    bPoint = aPoint;
    G4ThreeVector bPosition(0.,0.,theStep);
    bPoint->SetPosition(bPosition);
    step->SetPostStepPoint(bPoint);
    step->SetStepLength(theStep);

  // --------- Test the DoIt 
    G4int nElectrons = 0;
    G4int nPositrons = 0;
    G4int nPhotons = 0;
    G4double rmax = 0.0;
    G4double de = 0.0;
    G4double de2 = 0.0;

    G4cout << "dProcess= " << dProcess << "  cdProcess= " << cdProcess << G4endl;

    timer = new G4Timer();
    timer->Start();

    for (G4int iter=0; iter<nEvt; iter++) {

      gTrack->SetStep(step); 
 
      if(verbose) {
        G4cout  <<  "Iteration = "  <<  iter 
	        << "  -  Step Length = " 
	        << step->GetStepLength()/mm << " mm "
	        << G4endl;
      }

      G4VParticleChange* dummy = 0;
      if(postDo) {
        if(dProcess)  dummy = dProcess->PostStepDoIt(*gTrack, *step);
        if(cdProcess) dummy = cdProcess->PostStepDoIt(*gTrack, *step);
      } else {
        dummy = cdProcess->AlongStepDoIt(*gTrack, *step);
      }
      G4ParticleChange* particleChange = (G4ParticleChange*) dummy;
      
      // Primary physical quantities 

      G4double energyChange = particleChange->GetEnergyChange();
      G4double deltaE = gEnergy - energyChange ;
      G4double dedx = deltaE / (step->GetStepLength());
      
      G4ThreeVector change = particleChange->CalcMomentum(energyChange,
		            *(particleChange->GetMomentumChange()),
			     part->GetPDGMass());

      G4double pxChange  = change.x();
      G4double pyChange  = change.y();
      G4double pzChange  = change.z();

      if(verbose) {
        G4cout << "---- Primary after the step ---- " << G4endl;
        G4cout << "---- Energy: " << energyChange/MeV << " MeV,  " 
	       << "(px,py,pz): ("
	       << pxChange/MeV << ","
	       << pyChange/MeV << "," 
	       << pzChange/MeV << ") MeV"
	       << G4endl;
        G4cout << "---- Energy loss (dE) = " << deltaE/keV << " keV;" 
               << "Stopping power (dE/dx)=" << dedx*mm/keV << " keV/mm" 
               << "; rmax(mm)= " << rmax << G4endl;
      }

      // Primary
 
      G4int nsec = particleChange->GetNumberOfSecondaries();
      /*
      if(ntuple1) {
        ntuple1->column("epri", gEnergy/MeV);
        ntuple1->column("efin", energyChange/MeV);
        ntuple1->column("dedx", dedx*mm/MeV);
        ntuple1->column("nsec", nsec);
        ntuple1->column("nele", nElectrons);
        ntuple1->column("npho", nPhotons);
        ntuple1->dumpData(); 
      }
      */
      de  += deltaE;
      de2 += deltaE*deltaE; 
      
      // Secondaries physical quantities 
      
      if(usepaw) {
        hist[2]->fill((float)nsec, 1.0);
        hist[3]->fill(particleChange->GetLocalEnergyDeposit()/MeV, 1.0);
      }      

      for (G4int i = 0; i<nsec; i++) {
	  // The following two items should be filled per event, not
	  // per secondary; filled here just for convenience, to avoid
	  // complicated logic to dump ntuple when there are no secondaries
	  
        G4Track* finalParticle = particleChange->GetSecondary(i) ;
	  
        G4double e    = finalParticle->GetTotalEnergy();
        G4double eKin = finalParticle->GetKineticEnergy();
        G4double px   = (finalParticle->GetMomentum()).x();
        G4double py   = (finalParticle->GetMomentum()).y();
        G4double pz   = (finalParticle->GetMomentum()).z();
        //G4double theta= (finalParticle->GetMomentum()).theta();
        G4double p    = sqrt(px*px + py*py + pz*pz);

        if (eKin > gEnergy) {
	    G4cout << "WARNING: eFinal > eInit in event #" << iter << G4endl;
        }

        G4String partName = finalParticle->GetDefinition()->GetParticleName();
        if(verbose) {
	  G4cout  << "==== Final " 
		  <<  partName  <<  " "  
		  << "E= " <<  e/MeV  <<  " MeV,  " 
		  << "eKin: " <<  eKin/MeV  <<  " MeV, " 
		  << "(px,py,pz): ("
		  <<  px/MeV  <<  "," 
		  <<  py/MeV  <<  ","
		  <<  pz/MeV  << ") MeV," 
                  << " p= " << p << " MeV" 
		  <<  G4endl;   
	}
	  
        if(usepaw) {
	  hist[0]->fill(eKin/gEnergy, 1.0);
          hist[1]->fill(p/MeV, 1.0);
	}
	  
        G4int partType = 0;
        if (partName == "e-") {
	  partType = 1;
          nElectrons++;
	   
	} else if (partName == "e+") { 
	  partType = 2;
          nPositrons++;
	   
	} else if (partName == "gamma") { 
	  partType = 0;
          nPhotons++;
        }
	
	// Fill the secondaries ntuple
	/*
        if(ntuple2) {
          ntuple2->column("eprimary",gEnergy);
          ntuple2->column("px", px);
          ntuple2->column("py", py);
          ntuple2->column("pz", pz);
          ntuple2->column("p", p);
          ntuple2->column("e", e);
          ntuple2->column("theta", theta);
          ntuple2->column("ekin", eKin);
          ntuple2->column("type", partType);  
	  ntuple2->dumpData(); 
	}
	*/
	delete particleChange->GetSecondary(i);
      }
	          
      particleChange->Clear();      
    
    }
    G4cout << "###### Statistics:" << G4endl;
    G4cout << "Average number of secondary electrons= " 
           << (G4double)nElectrons/(G4double)nEvt << G4endl;
    G4cout << "Average number of secondary positrons= " 
             << (G4double)nPositrons/(G4double)nEvt << G4endl;
    G4cout << "Average number of secondary photons= " 
           << (G4double)nPhotons/(G4double)nEvt << G4endl;
    G4double x = de/(G4double)nEvt;
    G4double y = de2/(G4double)nEvt - x*x;
    if(0.0 < y) y = sqrt(y); 
    G4cout << "Average energy deposition(MeV)= " 
           << x/MeV << " +- " << y/MeV << G4endl;

    timer->Stop();
    G4cout << "  "  << *timer << G4endl;
    delete timer;
  
    if(usepaw) {
      tree->commit();
      G4std::cout << "Closing the tree..." << G4std::endl;
      tree->close();
      G4cout << "# hbook is writed" << G4endl;
    }

    G4cout << "###### End of run # " << run << "     ######" << G4endl;
    
  } while(end);
  G4cout << "###### End of test #####" << G4endl;
}












