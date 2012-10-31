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
//
#include <iostream>
#include <iomanip>
#include <fstream>

#include "BeamTestSiliconMonitorHit.hh"
#include "G4SystemOfUnits.hh"
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
