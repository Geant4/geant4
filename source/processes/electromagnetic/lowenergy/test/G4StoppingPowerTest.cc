// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      File name:     G4StoppingPowerTest
//                     This test provide cross sections 
//                     tests for electromagnetic processes. The input
//                     data have to be describe in ASCII file
//
//      Author:        V.Ivanchenko 
// 
//      Creation date: 23 May 2001
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
#include "G4VParticleChange.hh"

#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyRayleigh.hh"
#include "G4hLowEnergyIonisation.hh"

#include "G4VeEnergyLoss.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eplusAnnihilation.hh"

//#include "G4VeLowEnergyLoss.hh"

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
#include "G4ForceCondition.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4GRSVolume.hh"

#include "G4UnitsTable.hh"
#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"

HepTupleManager* hbookManager;

int main(int argc,char** argv)
{

  // -------------------------------------------------------------------
  // Setup

  G4int  nPart       =-1;
  G4String  nameMat  = "Si";
  G4int  nProcess    = 3;
  G4bool usepaw      = false;
  G4bool lowE        = true;
  G4int verbose      = 0;
  G4double emin      = 0.01*MeV;
  G4double emax      = 100.0*MeV;
  G4int nbin         = 1000;
  G4String hFile     = "";
  G4double theStep   = 0.01*micrometer;
  G4double range     = 1.0*micrometer;
  G4double cutG      = 10.0*mm;
  G4double cutE      = 10.0*mm;
  G4Material* material = 0; 
  G4String name[3] = {"Ionisation", "Bremsstrahlung",
                      "Ionisation+Bremsstrahlung"};


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
 
  G4Material* ma[14];
  ma[0] = new G4Material("Be",    4.,  9.01*g/mole, 1.848*g/cm3);
  ma[1] = new G4Material("Graphite",6., 12.00*g/mole, 2.265*g/cm3 );
  ma[2] = new G4Material("Al", 13., 26.98*g/mole, 2.7 *g/cm3);
  ma[3] = new G4Material("Si",   14., 28.055*g/mole, 2.33*g/cm3);
  /*
  ma[4] = new G4Material("LAr",   18., 39.95*g/mole, 1.393*g/cm3);
  ma[5] = new G4Material("Fe",      26., 55.85*g/mole, 7.87*g/cm3);
  ma[6] = new G4Material("Cu",    29., 63.55*g/mole, 8.96*g/cm3);
  ma[7] = new G4Material("W", 74., 183.85*g/mole, 19.30*g/cm3);
  ma[8] = new G4Material("Pb",82., 207.19*g/mole, 11.35*g/cm3);
  ma[9] = new G4Material("U", 92., 238.03*g/mole, 18.95*g/cm3);

  G4Element*   H  = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   O  = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
  G4Element*   C  = new G4Element ("Carbon"  , "C", 6. , 12.00*g/mole);
  G4Element*  Cs  = new G4Element ("Cesium"  , "Cs", 55. , 132.905*g/mole);
  G4Element*   I  = new G4Element ("Iodide"  , "I", 53. , 126.9044*g/mole);

  ma[10] = new G4Material("O2", 8., 16.00*g/mole, 1.1*g/cm3);
  ma[11] = new G4Material ("Water" , 1.*g/cm3, 2);
  ma[11]->AddElement(H,2);
  ma[11]->AddElement(O,1);

  ma[12] = new G4Material ("Ethane" , 0.4241*g/cm3, 2);
  ma[12]->AddElement(H,6);
  ma[12]->AddElement(C,2);
  
  ma[13] = new G4Material ("CsI" , 4.53*g/cm3, 2);
  ma[13]->AddElement(Cs,1);
  ma[13]->AddElement(I,1);
  */ 
  static const G4MaterialTable* theMaterialTable = 
               G4Material::GetMaterialTable();

  G4int nMaterials = theMaterialTable->length();
  G4cout << "Available materials are: " << G4endl;
  G4int mat;
  for (mat = 0; mat < nMaterials; mat++) {
    G4cout << mat << ") " << ma[mat]->GetName() << G4endl;
  }

  G4cout << "Available processes are: " << G4endl;
  for (mat = 0; mat < 2; mat++) {
    G4cout << mat << ") " << name[mat] << G4endl;
  }

  // Particle definitions

  G4ParticleDefinition* gamma = G4Gamma::GammaDefinition();
  G4ParticleDefinition* electron = G4Electron::ElectronDefinition();
  G4ParticleDefinition* positron = G4Positron::PositronDefinition();
  G4ParticleDefinition* proton   = G4Proton::ProtonDefinition();
  G4ParticleDefinition* antiproton = G4AntiProton::AntiProtonDefinition();
  G4ParticleDefinition* part = electron;

  G4hLowEnergyIonisation* hionle = 0;
  G4LowEnergyIonisation* eionle = 0;
  G4LowEnergyBremsstrahlung* ebrle = 0;
  G4hIonisation* hionst = 0;
  G4eIonisation* eionst = 0;
  G4eBremsstrahlung* ebrst = 0;

  G4cout  <<  "Process is initialized"  <<  G4endl;; 

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
  G4cout << "#emin(MeV)" << G4endl;
  G4cout << "#emax(MeV)" << G4endl;
  G4cout << "#nbin" << G4endl;
  G4cout << "#cutG(mm)" << G4endl;
  G4cout << "#cutE(mm)" << G4endl;
  G4cout << "#range(mm)" << G4endl;
  G4cout << "#step(mm)" << G4endl;
  G4cout << "#material" << G4endl;
  G4cout << "#process" << G4endl;
  G4cout << "#domain" << G4endl;
  G4cout << "#paw" << G4endl;
  G4cout << "#verbose" << G4endl;
  G4cout << "#run" << G4endl;
  G4cout << "#exit" << G4endl;

  G4ProcessManager *elecManager, *positManager, 
                   *protManager, *aprotManager;


  string line, line1;
  G4bool end = true;
  for(G4int run=0; run<100; run++) {
    do {
      (*fin) >> line;
      G4cout << "Next line " << line << G4endl;
      if(line == "#particle") {
        (*fin) >> nPart;
      } else if(line == "#emin(MeV)") {
        (*fin) >> emin;
        emin *= MeV;
      } else if(line == "#emax(MeV)") {
        (*fin) >> emax;
        emax *= MeV;
      } else if(line == "#nbin") {
        (*fin) >> nbin;
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
    if(nProcess < 0 || nProcess > 2) {
      G4cout << "Process #" << nProcess 
             << " is absent in the list of processes: Exit" << G4endl;
      break;
    }


    for (mat = 0; mat < nMaterials; mat++) {
      material = ma[mat];
      if(nameMat == material->GetName()) break;
    }

    G4cout << "The particle: " << part->GetParticleName() << G4endl;
    G4cout << "The material: " << material->GetName() << G4endl;
    G4cout << "The cut on e-:" << cutE/mm << " mm" << G4endl;
    G4cout << "The cut on g: " << cutG/mm << " mm" << G4endl;
    G4cout << "The step:     " << theStep/mm << " mm" << G4endl;

    // -------------------------------------------------------------------
    // ---- HBOOK initialization

    hbookManager = new HBookFile(hFile, 58);
    assert (hbookManager != 0);
  
    // ---- Book a histogram and ntuples
    G4cout << "Hbook file name: <" 
           << ((HBookFile*) hbookManager)->filename() << ">" << G4endl;
    G4double emin10 = log10(emin/MeV);
    G4double emax10 = log10(emax/MeV);
    G4double bin = (emax10 - emin10) / (G4double)nbin;

    HepHistogram* hist[2];
    hist[0] = hbookManager->histogram("Stopping power (MeV*cm**2/g)", 
                                     nbin,emin10,emax10);
    hist[1] = hbookManager->histogram("Step limit (mm)", 
                                     nbin,emin10,emax10);
    //    assert (hDebug != 0);  
    G4cout<< "Histograms is initialised" << G4endl;

    gamma->SetCuts(cutG);
    electron->SetCuts(cutE);
    positron->SetCuts(cutE);
  
    // Processes - all new 
    G4bool success = false;


    if(lowE) {

      elecManager = new G4ProcessManager(electron);
      electron->SetProcessManager(elecManager);
      eionle = new G4LowEnergyIonisation();
      eionle->SetEnlossFluc(false);
      ebrle = new G4LowEnergyBremsstrahlung();
      ebrle->SetEnlossFluc(false);
      elecManager->AddProcess(eionle);
      elecManager->AddProcess(ebrle);
      eionle->BuildPhysicsTable(*electron);
      ebrle->BuildPhysicsTable(*electron);
      if(nPart == 1) {
        if(nProcess == 2) {
          success = true;
	}

      } else if (nPart == 3) {
        protManager = new G4ProcessManager(proton);
        proton->SetProcessManager(protManager);
        hionle = new G4hLowEnergyIonisation();
        hionle->SetEnlossFluc(false);
        hionle->SetNuclearStoppingOff();
        protManager->AddProcess(hionle);
        hionle->BuildPhysicsTable(*proton);
        success = true;

      } else if (nPart == 4) {
        aprotManager = new G4ProcessManager(antiproton);
        antiproton->SetProcessManager(aprotManager);
        hionle = new G4hLowEnergyIonisation();
        hionle->SetEnlossFluc(false);
        hionle->SetNuclearStoppingOff();
        aprotManager->AddProcess(hionle);
        hionle->BuildPhysicsTable(*antiproton);
        success = true;
      }

    } else {

      elecManager = new G4ProcessManager(electron);
      electron->SetProcessManager(elecManager);
      eionst = new G4eIonisation();
      eionst->SetEnlossFluc(false);
      elecManager->AddProcess(eionst);
      ebrst = new G4eBremsstrahlung();
      ebrst->SetEnlossFluc(false);
      elecManager->AddProcess(ebrst);
      eionst->BuildPhysicsTable(*electron);
      ebrst->BuildPhysicsTable(*electron);
      if(nPart == 1) {
        if(nProcess == 2) {
          success = true;
	}
      } else if(nPart == 2) {
        positManager = new G4ProcessManager(positron);
        positron->SetProcessManager(positManager);
        if(nProcess == 0) {
          eionst = new G4eIonisation();
          eionst->SetEnlossFluc(false);
          positManager->AddProcess(eionst);
          eionst->BuildPhysicsTable(*positron);
          success = true;
        } else if(nProcess == 1) {
          ebrst = new G4eBremsstrahlung();
          ebrst->SetEnlossFluc(false);
          positManager->AddProcess(ebrst);
          ebrst->BuildPhysicsTable(*positron);
          success = true;
        } else if(nProcess == 2) {
          eionst = new G4eIonisation();
          ebrst = new G4eBremsstrahlung();
          eionst->SetEnlossFluc(false);
          ebrst->SetEnlossFluc(false);
          positManager->AddProcess(eionst);
          positManager->AddProcess(ebrst);
          eionst->BuildPhysicsTable(*positron);
          ebrst->BuildPhysicsTable(*positron);
          success = true;
	}

      } else if (nPart == 3) {
        protManager = new G4ProcessManager(proton);
        proton->SetProcessManager(protManager);
        hionst = new G4hIonisation();
        hionst->SetEnlossFluc(false);
        protManager->AddProcess(hionst);
        hionst->BuildPhysicsTable(*proton);
        success = true;

      } else if (nPart == 4) {
        aprotManager = new G4ProcessManager(antiproton);
        antiproton->SetProcessManager(aprotManager);
        hionst = new G4hIonisation();
        hionst->SetEnlossFluc(false);
        aprotManager->AddProcess(hionst);
        hionst->BuildPhysicsTable(*antiproton);
        success = true;
      }
    }

    if(success) G4cout  <<  "Physics tables are built"  <<  G4endl;  
    else        G4cout  <<  "Physics tables are not built!!!"  <<  G4endl;  


    // Create a DynamicParticle  
  
    G4ParticleMomentum gDir(initX,initY,initZ);
    G4double gEnergy = emax;
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

    for (G4int iter=0; iter<nbin; iter++) {

      G4double le = emin10 + ((G4double)iter + 0.5)*bin;
      G4double  e = pow(10.0,le) * MeV;
      gTrack->SetStep(step); 
      gTrack->SetKineticEnergy(e); 
      G4double x = 0.0;
      G4VParticleChange* aChange = 0;

      if(lowE) {

        if(nPart == 1) {
          if(nProcess == 0) {
            x = eionle->GetContinuousStepLimit(*gTrack,theStep,theStep,safety);
            aChange = eionle->AlongStepDoIt(*gTrack,*step);

          } else if(nProcess == 1) {
            x = ebrle->GetContinuousStepLimit(*gTrack,theStep,theStep,safety);
            aChange = ebrle->AlongStepDoIt(*gTrack,*step);
          } else if(nProcess == 2) {
            x = eionle->GetContinuousStepLimit(*gTrack,theStep,theStep,safety);
            aChange = eionle->AlongStepDoIt(*gTrack,*step);
 	  }

        } else if (nPart == 3) {
          x = hionle->GetContinuousStepLimit(*gTrack,theStep,theStep,safety);
          aChange = hionle->AlongStepDoIt(*gTrack,*step);

        } else if (nPart == 4) {
          x = hionle->GetContinuousStepLimit(*gTrack,theStep,theStep,safety);
          aChange = hionle->AlongStepDoIt(*gTrack,*step);
	}

      } else {

        if(nPart == 1) {
          if(nProcess == 0) {
            x = eionst->GetContinuousStepLimit(*gTrack,theStep,theStep,safety);
            aChange = eionst->AlongStepDoIt(*gTrack,*step);

          } else if(nProcess == 1) {
            x = ebrst->GetContinuousStepLimit(*gTrack,theStep,theStep,safety);
            aChange = ebrst->AlongStepDoIt(*gTrack,*step);
          } else if(nProcess == 2) {
            x = eionst->GetContinuousStepLimit(*gTrack,theStep,theStep,safety);
            aChange = eionst->AlongStepDoIt(*gTrack,*step);
 	  }

        } else if(nPart == 2) {
          if(nProcess == 0) {
            x = eionst->GetContinuousStepLimit(*gTrack,theStep,theStep,safety);
            aChange = eionst->AlongStepDoIt(*gTrack,*step);

          } else if(nProcess == 1) {
            x = ebrst->GetContinuousStepLimit(*gTrack,theStep,theStep,safety);
            aChange = ebrst->AlongStepDoIt(*gTrack,*step);
          } else if(nProcess == 2) {
            x = eionst->GetContinuousStepLimit(*gTrack,theStep,theStep,safety);
            aChange = eionst->AlongStepDoIt(*gTrack,*step);
 	  }

        } else if (nPart == 3) {
          x = hionst->GetContinuousStepLimit(*gTrack,theStep,theStep,safety);
          aChange = hionst->AlongStepDoIt(*gTrack,*step);

        } else if (nPart == 4) {
          x = hionst->GetContinuousStepLimit(*gTrack,theStep,theStep,safety);
          aChange = hionst->AlongStepDoIt(*gTrack,*step);
        }
      }

      G4double de = aChange->GetLocalEnergyDeposit();
      G4int n = aChange->GetNumberOfSecondaries();

      //G4cout << " de(MeV) = " << de/MeV << " n= " << n << G4endl;


      if(n > 0) {
        for(G4int i=0; i<n; i++) {
          de += (aChange->GetSecondary(i))->GetKineticEnergy();
          G4cout << "add " 
                 << ((aChange->GetSecondary(i))->GetKineticEnergy())/eV
                 << " eV" << G4endl;
	}
      }
      G4double st = de/(theStep*(material->GetDensity()));
      st *= gram/(cm*cm*MeV); 
 
      if(verbose) {
          G4cout  <<  "Iteration = "  <<  iter 
	          << "  E = " << e/MeV << " MeV; StepLimit= "
	          << x/mm << " mm; de= " 
                  << de/eV << " eV; dE/dx= "
                  << st << " MeV*cm^2/g" <<  G4endl;
      }

      if(x > 1000.0*meter) x = 1000.0*meter;      

      hist[0]->accumulate(le,st);
      hist[1]->accumulate(le,x/mm);
    }
    if(usepaw)hbookManager->write();
    G4cout << "# hbook is writed" << G4endl;
    delete hbookManager;    
    G4cout << "# hbook is deleted" << G4endl;
    G4cout << "###### End of run # " << run << "     ######" << G4endl;
    

  } while(end);
  G4cout << "###### End of test #####" << G4endl;
}












