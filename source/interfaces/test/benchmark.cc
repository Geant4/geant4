// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: benchmark.cc,v 1.1 1999-01-08 16:32:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------------
// 
// Test routine for combining "events+tracks", "tracking","geometry"
// and "particle+matter"
//
// Calorimeter geometry
//
//   Arranged by R.Kiuchi 09. Oct. 1995
//
//   Original source written by Makoto Asai   5/Sep/95
//   revised for New G4ParticleDefinition by H.Kurashige  21/Apr/96
//   revised for New G4ParticleGun by H.Kurashige  3/May/96
//   revised for interactive version by M.Asai 11/Jun/96
//   revised for tag 0.5 by G.Cosmo  12/Feb/97
//   revised for alpha01-tag by G.Cosmo  8/Apr/97
// 

#include "G4ios.hh"
#include <stdlib.h>

#include "G4ThreeVector.hh"

#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"

#include "G4GeometryManager.hh"
#include "G4TransportationManager.hh"

#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4ProcessManager.hh"
#include "G4MuDeltaRay.hh"
#include "G4MuPairProduction.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuEnergyLoss.hh"
#include "G4Transportation.hh"
#include "G4Geantino.hh"
#include "G4Electron.hh"
#include "G4MuonPlus.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"

G4VPhysicalVolume* BuildCalorimeter(G4Material* Air,
				    G4Material* Lead,
				    G4Material* Alluminium)
{
    G4double offset=22.5*cm, xTlate, yTlate;
    G4int i,j,copyNo;

    G4Box *myWorldBox= new G4Box ("WBox",10000*cm, 10000*cm, 10000*cm);
    G4Box *myCalBox = new G4Box ("CBox",1500*cm, 1500*cm, 1000*cm);
    G4Tubs *myTargetTube = new G4Tubs ("TTube",0*cm, 22.5*cm, 1000*cm, 0.*deg, 360.*deg);

    G4LogicalVolume *myWorldLog=new G4LogicalVolume(myWorldBox,Air,
						    "WLog", 0, 0, 0);
    G4LogicalVolume *myCalLog=new G4LogicalVolume(myCalBox,Alluminium,
						  "CLog", 0, 0, 0);
    G4LogicalVolume *myTargetLog=new G4LogicalVolume(myTargetTube,Lead,
						     "TLog", 0, 0, 0);

    G4PVPlacement *myWorldPhys=new G4PVPlacement(0,G4ThreeVector(),
						 "WPhys",
						 myWorldLog,
						 0,false,0);
    G4PVPlacement *myCalPhys=new G4PVPlacement(0,G4ThreeVector(),
					       "CalPhys",
					       myCalLog,
					       myWorldPhys,false,0);

    G4String tName1("TPhys1");	// Allow all target physicals to share
				// same name (delayed copy)
    copyNo=0;
    for (j=1;j<=25;j++)
      {
	yTlate = -1000.0*cm - 40.0*cm + j*80.0*cm;
	
	for (i=1;i<=50;i++)
	  {
	    copyNo++;
	    xTlate = -1000.0*cm - 20.0*cm + i*45.0*cm - offset;
	    G4PVPlacement *myTargetPhys=new G4PVPlacement(0,G4ThreeVector(xTlate,yTlate,0*cm),
							  tName1,
							  myTargetLog,
							  myCalPhys,
							  false,
							  copyNo);
	  }
      }

    G4String tName2("TPhys2");	// Allow all target physicals to share
				// same name (delayed copy)
    copyNo=0;
    for (j=1;j<=26;j++)
	{
	  yTlate = -1000.0*cm - 80.0*cm + j*80.0*cm;
	  for (i=1;i<=50;i++)
	    {
	      copyNo++;
	      xTlate = -1000.0*cm - 20.0*cm + i*45.0*cm;
	      G4PVPlacement *myTargetPhys=new G4PVPlacement(0,G4ThreeVector(xTlate,yTlate,0*cm),
							    tName2,
							    myTargetLog,
							    myCalPhys,
							    false,
							    copyNo);
	    }
	}

    return myWorldPhys;
}

int main()
{
//-------- set output format-------
    G4cout.setf( ios:: scientific, ios :: floatfield );

//--------- Material definition ---------

  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxigen", symbol="O", iz=8., a);

  a = 26.98*g/mole;
  density = 2.7*g/cm3;
  G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);
  a = 55.85*g/mole;
  density = 7.87*g/cm3;
  G4Material* Fe = new G4Material(name="Iron", z=26., a, density);
  a = 207.19*g/mole;
  density = 11.35*g/cm3;
  G4Material* Pb = new G4Material(name="Lead", z=82., a, density);
  density = 1.29e-03*g/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);

// Geometry definitions

    G4VPhysicalVolume *myTopNode;

    G4cout << "Calorimeter volume Performance Test" << endl;
    myTopNode=BuildCalorimeter(Air,Pb,Al);	// Build the geometry
    G4GeometryManager::GetInstance()->CloseGeometry();

    G4TransportationManager::GetTransportationManager()
      ->GetNavigatorForTracking()->SetWorldVolume( myTopNode );

// Particle and Process definitions

// Particle definitions

    G4ParticleDefinition* Gamma    = G4Gamma::Gamma();
    G4ParticleDefinition* Geantino = G4Geantino::Geantino();
    G4ParticleDefinition* MuonPlus = G4MuonPlus::MuonPlus();
    G4ParticleDefinition* Electron = G4Electron::Electron();
    G4ParticleDefinition* Positron = G4Positron::Positron();
    G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();

// Process definitions

   G4ProcessManager* pmanager;
   G4Transportation* theTransportationProcess= new G4Transportation();

   G4ParticleTable::G4PTblDicIterator *piter = theParticleTable->GetIterator();
   piter -> reset();
   while( (*piter)() ){
     pmanager = piter->value()->GetProcessManager();
     // add Transportation
     pmanager -> AddContinuousProcess(theTransportationProcess, 1);
   }

/* !!!!!!!!!!!!!!!!!!!  Muons temporary commented out !!!!!!!!!!!!!!!!!
    G4MuDeltaRay theDeltaRayProcess;
    G4MuPairProduction thePairProductionProcess;
    G4MuBremsstrahlung theBremsstrahlungProcess;

    G4MuEnergyLoss theEnergyLossProcess;

    G4ProcessManager* muonProcessManager= MuonPlus->GetProcessManager();
    G4ProcessManager* geantinoProcessManager= Geantino->GetProcessManager();
    G4ProcessManager* electronProcessManager= Electron->GetProcessManager();

    muonProcessManager->AddDiscreteProcess(&theDeltaRayProcess);
    muonProcessManager->AddDiscreteProcess(&thePairProductionProcess); 
    muonProcessManager->AddDiscreteProcess(&theBremsstrahlungProcess);

    muonProcessManager->AddContinuousProcess(&theEnergyLossProcess);

// Set the cut and Build cross-sections table for physics processes

    Gamma->SetCuts(1.*cm);
    Positron->SetCuts(1.*cm);
    Electron->SetCuts(1.*cm);
    MuonPlus->SetCuts(1.*cm);
 !!!!!!!!!!!!!!!!!!!  Muons temporary commented out !!!!!!!!!!!!!!!!! */

// Event manager

    G4EventManager * eventManager = new G4EventManager;

// Event generator

    G4int numPrimaries =1;    //1 particle per event
    G4ParticleGun * particleGun = new G4ParticleGun(numPrimaries);

// Set navigator from transportation manager

    eventManager->
      set_Navigator(G4TransportationManager::GetTransportationManager()->
                    GetNavigatorForTracking() );

// Prepare UI session

    G4UImanager * UI = G4UImanager::GetUIpointer();
    G4UIterminal * terminal = new G4UIterminal;
    UI->SetSession( terminal );

// OK. Let's start benchmarks............

    UI->Interact("G4> ");

// End of run business

    G4GeometryManager::GetInstance()->OpenGeometry();

    delete Air; delete Pb; delete Fe; delete Al;
    delete elO; delete elN;
    delete theTransportationProcess;
    delete eventManager;
    delete particleGun;

    return EXIT_SUCCESS;
}
