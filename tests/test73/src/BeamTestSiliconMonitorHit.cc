//
#include "BeamTestSiliconMonitorHit.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include "G4UnitsTable.hh"

G4Allocator<BeamTestSiliconMonitorHit> BeamTestSiliconMonitorHitAllocator;

BeamTestSiliconMonitorHit::BeamTestSiliconMonitorHit()
{
  /* fIPD = 0;
   fIKEnergy = 0.0;
   fIPosition = G4ThreeVector( 0.0,0.0,0.0 );
   fIMomentumD = G4ThreeVector( 0.0,0.0,0.0 );*/
   fEPD = 0;
   fEKEnergy = 0.0;
   fEPosition = G4ThreeVector( 0.0,0.0,0.0 );
   fEMomentumD = G4ThreeVector( 0.0,0.0,0.0 );
   fEMomentum = G4ThreeVector( 0.0,0.0,0.0 );
    fChamberNum = 0;
    fTrackId = -1;
    fEnergy = 0;
    fStatus = fUndefined;
    fStepLength = 0;
}



BeamTestSiliconMonitorHit::~BeamTestSiliconMonitorHit() {}



void BeamTestSiliconMonitorHit::Draw() {}



void BeamTestSiliconMonitorHit::Print()
{
	/*
	   G4cout << "Incidence Particle Name and Kinetic Energy " << fIPD->GetParticleName() 
	   << " " << fIKEnergy/MeV << " MeV" << G4endl; 
	   G4cout << "Insidence position in silicon monitor " << fIPosition/mm << " in mm" << G4endl; 
	   G4cout << "Incidence Direction " << fIMomentumD << G4endl;
	 */
	G4cout << "Exiting Particle Name and Kinetic Energy " << fEPD->GetParticleName() 
		<< " " << fEKEnergy/MeV << " MeV" << G4endl; 
	G4cout << "Exiting position in silicon monitor " << fEPosition/mm << " in mm" << G4endl; 
	 
	//std::cout << "Exiting X position in silicon monitor " << std::setprecision(5) << fEPosition.x()/mm << " in mm" << std::endl; 
	G4cout << "Exiting Momentum Direction " << fEMomentumD << G4endl;
	G4cout << "Exiting Position " << fEMomentum << G4endl;
    G4cout << "Chamber Number "<<fChamberNum<<G4endl;
    G4cout << "Track ID "<<fTrackId<<G4endl;
	//std::fstream out;
   	//out.open("test.csv", std::ios::out | std::ios::app);
	//out << fEPosition.x()/mm << ", " << fEPosition.y()/mm << ", " << fEPosition.z()/mm << std::endl;
	//std::ofstream outfile ("test.cvs");
	//outfile << std::setprecision(10) << fEPosition.x()/mm;
	//outfile.close();
	
}
