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
// 
#include <CLHEP/Units/SystemOfUnits.h>

#include "BeamTestPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4RunManager.hh"
#include "BeamTestRunAction.hh"
#include "BeamTestEventAction.hh"

BeamTestPrimaryGeneratorAction::BeamTestPrimaryGeneratorAction(/*Parameters* parameter*/)
{
    //Set defaults
	G4int n_particle = 1;
	particleGun  = new G4ParticleGun(n_particle);

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;

	G4ParticleDefinition* particle = particleTable->FindParticle(particleName="mu-");

	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
	//particleGun->SetParticleEnergy(10.*GeV);
}

BeamTestPrimaryGeneratorAction::~BeamTestPrimaryGeneratorAction()
{
	delete particleGun;
}

void BeamTestPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    /*
    BeamTestConversion A(parameter);
	G4cout << G4endl;
	G4cout << "Angle " << xzAngle(A) << " " << A.getP() << " " << A.getPT() << " " << "New Pz " << pZ(A) << G4endl;
	G4cout << G4endl;
    
	//G4double pT = A.getPT();
	G4double p = A.getP();
	G4double theta = xzAngle(A);
	//G4double pz = pZ(A);
    
	//this->Position(parameter);
	//particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
	particleGun->SetParticleMomentumDirection(G4ThreeVector(std::tan(theta),0.,1.0));
	//particleGun->SetParticleEnergy(10.*GeV);
	particleGun->SetParticleEnergy(p*GeV);
	// This forces the particle gun to be exactly 10cm from the first detector plane
	// We can also set the number of Particles
	//particleGun->SetNumberOfParticles(100);
     */
    
    //Get info on kinematics from particleGun object: now updated to 
    // values found in macro file
    if ( anEvent->GetEventID() == 0 ) { //First event of the run, setup stuff...
        //Get momentum direction
        G4ThreeVector momdir = particleGun->GetParticleMomentumDirection();
        //Position of primary depends on P
        G4double dx = 6.*CLHEP::mm;
        G4double dz = dx/std::tan(momdir.theta());
        //G4cout << "pT " << pT << "  and tan(theta_xz) = " << tan(xzAngle(A)) << G4endl;
        G4ThreeVector Position(-dx,0.,-dz);
        particleGun->SetParticlePosition(Position);
    
        //Configure other user actions that require information on beam parameters
        G4UserRunAction* theRunAction = const_cast< G4UserRunAction* > (G4RunManager::GetRunManager()->GetUserRunAction());
        //Get mometum magnitude
        G4double myp = particleGun->GetParticleMomentum();
        //Pt.. Remember only along x... 
        G4double mypt = myp * std::sin(momdir.theta());
        static_cast<BeamTestRunAction*>(theRunAction)->SetpT(mypt);
        static_cast<BeamTestRunAction*>(theRunAction)->Initialize();
        G4UserEventAction* theEventAction = const_cast<G4UserEventAction*>(G4RunManager::GetRunManager()->GetUserEventAction());
        G4ThreeVector mymom = momdir;
        mymom.setMag( myp );
        static_cast<BeamTestEventAction*>(theEventAction)->Initialize( mymom, Position );
        G4cout<<"Primaries "<<particleGun->GetParticleDefinition()->GetParticleName()<<" will be generated with Momentum: "<<G4BestUnit(mymom,"Energy")<<"/c , at Postion: "<<G4BestUnit(Position,"Length")<<G4endl;
    }
    //Generate a primary (i.e. add it to the stack of particles)
    particleGun->GeneratePrimaryVertex(anEvent);

}

/*void BeamTestPrimaryGeneratorAction::Position(Parameters* parameter)
{
	
}*/
