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
// $Id: G4ComplexTest.cc,v 1.7 2001-10-08 16:36:07 vnivanch Exp $
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
#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"

//typedef G4LowEnergyBremsstrahlungIV G4LowEnergyBremsstrahlung;
//typedef G4eLowEnergyIonisationIV G4LowEnergyIonisation;

int main(int argc,char** argv)
{

  // -------------------------------------------------------------------
  // Setup

  HepTupleManager* hbookManager;
  HepTuple* ntuple1 = 0; 
  HepTuple* ntuple2 = 0; 
  HepHistogram* h[4] = {0,0,0,0};


  G4int  nEvt        = 100;
  G4int  nPart       =-1;
  G4String  nameMat  = "Si";
  G4int  nProcess    = 3;
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

  G4Material* Be = new G4Material("Be",    4.,  9.01*g/mole, 1.848*g/cm3);
  G4Material* Graphite = new G4Material("Graphite",6., 12.00*g/mole, 2.265*g/cm3 );
  G4Material* Al  = new G4Material("Al", 13., 26.98*g/mole, 2.7 *g/cm3);
  G4Material* Si  = new G4Material("Si",   14., 28.055*g/mole, 2.33*g/cm3);
  G4Material* LAr = new G4Material("LAr",   18., 39.95*g/mole, 1.393*g/cm3);
  G4Material* Fe  = new G4Material("Fe",      26., 55.85*g/mole, 7.87*g/cm3);
  G4Material* Cu  = new G4Material("Cu",    29., 63.55*g/mole, 8.96*g/cm3);
  G4Material*  W  = new G4Material("W", 74., 183.85*g/mole, 19.30*g/cm3);
  G4Material* Pb  = new G4Material("Pb",      82., 207.19*g/mole, 11.35*g/cm3);
  G4Material*  U  = new G4Material("U", 92., 238.03*g/mole, 18.95*g/cm3);

  G4Element*   H  = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   O  = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
  G4Element*   C  = new G4Element ("Carbon"  , "C", 6. , 12.00*g/mole);
  G4Element*  Cs  = new G4Element ("Cesium"  , "Cs", 55. , 132.905*g/mole);
  G4Element*   I  = new G4Element ("Iodide"  , "I", 53. , 126.9044*g/mole);

  G4Material*  maO = new G4Material("O2", 8., 16.00*g/mole, 1.1*g/cm3);

  G4Material* water = new G4Material ("Water" , 1.*g/cm3, 2);
  water->AddElement(H,2);
  water->AddElement(O,1);

  G4Material* ethane = new G4Material ("Ethane" , 0.4241*g/cm3, 2);
  ethane->AddElement(H,6);
  ethane->AddElement(C,2);
  
  G4Material* csi = new G4Material ("CsI" , 4.53*g/cm3, 2);
  csi->AddElement(Cs,1);
  csi->AddElement(I,1);

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  G4int nMaterials = theMaterialTable->length();

  G4cout << "Available materials are: " << G4endl;
  G4int mat;
  for (mat = 0; mat < nMaterials; mat++) {
    G4cout << mat << ") " << (*theMaterialTable)(mat)->GetName() << G4endl;
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
      material = (*theMaterialTable)(mat);
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

    hbookManager = new HBookFile(hFile, 58);
    //  assert (hbookManager != 0);
  
    // ---- Book a histogram and ntuples
    G4cout<< "Hbook file name: <" << ((HBookFile*) hbookManager)->filename() 
          << ">" << G4endl;

    /*
    // ---- primary ntuple ------
    ntuple1 = hbookManager->ntuple("Primary Ntuple");
  
    // ---- secondary ntuple ------
    ntuple2 = hbookManager->ntuple("Secondary Ntuple");
    */
    // ---- secondaries histos ----
    h[0] = hbookManager->histogram("Kinetic Energy (MeV)", 50,0.,gEnergy/MeV);
  
    h[1] = hbookManager->histogram("Momentum (MeV/c)", 50,0.,gEnergy*0.1/MeV);
  
    h[2] = hbookManager->histogram("Number of secondaries", 20,-0.5,19.5);
  
    h[3] = hbookManager->histogram("Energy deposition (MeV)", 50,0.,
                                     gEnergy*0.1/MeV);
    G4cout<< "Histograms is initialised" << G4endl;

    gamma->SetCuts(cutG);
    electron->SetCuts(cutE);
    positron->SetCuts(cutE);
  
    // Processes 

    G4VDiscreteProcess*            dProcess;
    G4VContinuousDiscreteProcess* cdProcess;
    dProcess = 0;
    cdProcess = 0;

    G4cout  <<  "Process is initialized"  <<  G4endl;; 

    if(lowE) {
      if(ionis) {
        elecLEion = new G4LowEnergyIonisation();
        elecLEbr = new G4LowEnergyBremsstrahlung();
        elecManager->AddProcess(elecLEion);
        elecManager->AddProcess(elecLEbr);
        elecLEbr->BuildPhysicsTable(*electron);
        elecLEion->BuildPhysicsTable(*electron);
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
        elecSTbr->BuildPhysicsTable(*electron);
        elecSTion->BuildPhysicsTable(*electron);
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
        G4GPILSelection* sel = 0;
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
      G4double pChange   = sqrt(pxChange*pxChange + pyChange*pyChange 
                                                  + pzChange*pzChange);
      
      G4double xChange = particleChange->GetPositionChange()->x();
      G4double yChange = particleChange->GetPositionChange()->y();
      G4double zChange = particleChange->GetPositionChange()->z();
      G4double thetaChange = particleChange->GetPositionChange()->theta();

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

      if(ntuple1) {
        ntuple1->column("epri", gEnergy/MeV);
        ntuple1->column("efin", energyChange/MeV);
        ntuple1->column("dedx", dedx*mm/MeV);
        ntuple1->column("nsec", nsec);
        ntuple1->column("nele", nElectrons);
        ntuple1->column("npho", nPhotons);
        ntuple1->dumpData(); 
      }
      de  += deltaE;
      de2 += deltaE*deltaE; 
      
      // Secondaries physical quantities 
      
      h[2]->accumulate((float)nsec, 1.0);
      h[3]->accumulate(particleChange->GetLocalEnergyDeposit()/MeV, 1.0);
      
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
        G4double theta= (finalParticle->GetMomentum()).theta();
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
	  
	h[0]->accumulate(eKin/MeV, 1.0);
        h[1]->accumulate(p/MeV, 1.0);
	  
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
  
    if(usepaw)hbookManager->write();
    G4cout << "# hbook is writed" << G4endl;
    delete hbookManager;    
    G4cout << "# hbook is deleted" << G4endl;
    G4cout << "###### End of run # " << run << "     ######" << G4endl;
    
  } while(end);
  G4cout << "###### End of test #####" << G4endl;
}












