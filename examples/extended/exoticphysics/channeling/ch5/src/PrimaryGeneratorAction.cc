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
// gpaterno, October 2025
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "Randomize.hh"

#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    //G4cout << "### PrimaryGeneratorAction instantiated ###" << G4endl;
    
    //instantiating the messenger
    fMessenger = new PrimaryGeneratorActionMessenger(this);
    
    //defining a Particle Gun
    G4int n_particle = 1;
    fGun = new G4ParticleGun(n_particle);

    //defining a GPS (used if fUseGPS==true)
    fGPS = new G4GeneralParticleSource(); 
    
    //default particle
    G4ParticleDefinition* particle = 
        G4ParticleTable::GetParticleTable()->FindParticle(fType);
    
    //set default values for GPS
    //Primary particles
    fGPS->SetParticleDefinition(particle);    
    //Position distribution
    G4SPSPosDistribution* vPosDist = 
        fGPS->GetCurrentSource()->GetPosDist();
    vPosDist->SetPosDisType("Plane");
    vPosDist->SetPosDisShape("Circle");
    vPosDist->SetRadius(0.5*CLHEP::mm);   
    vPosDist->SetCentreCoords(G4ThreeVector(0., 0., -5*CLHEP::cm));
    //Angular distribution
    G4SPSAngDistribution* vAngDist = 
        fGPS->GetCurrentSource()->GetAngDist();
    vAngDist->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));  
    //Energy distribution
    G4SPSEneDistribution* vEneDist = 
        fGPS->GetCurrentSource()->GetEneDist();
    vEneDist->SetEnergyDisType("Mono");
    vEneDist->SetMonoEnergy(fEnergy);
    
    //set default values for Particle Gun
    fGun->SetParticleDefinition(particle);
    fGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    fGun->SetParticleEnergy(fEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fGPS;
    delete fGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    if (fUseGPS) {
        fGPS->GeneratePrimaryVertex(anEvent); //override default through macro
    } else {              
        //particle position
        G4double x = G4RandGauss::shoot(fX, fSigmaX);
        G4double y = G4RandGauss::shoot(fY, fSigmaY);
        G4double z = G4RandGauss::shoot(fZ, fSigmaZ);
        fGun->SetParticlePosition(G4ThreeVector(x, y, z));

        //time, altrnative way to take into account the bunch length
        G4double t = fT;
        if (fSigmaT != 0 && fSigmaZ == 0) {
            t = G4RandGauss::shoot(fT, fSigmaT);
        }
        fGun->SetParticleTime(t);
        
        //particle Type
        G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
        fGun->SetParticleDefinition(particleTable->FindParticle(fType));
        
        //particle kinetic energy
        G4double K = G4RandGauss::shoot(fEnergy, fEnergy*fRelSigmaEnergy);
        fGun->SetParticleEnergy(K);   
           
        //particle direction      
        G4double xp = G4RandGauss::shoot(fXp, fSigmaXp); //xp=x'=px/pz=dx/dz
        G4double yp = G4RandGauss::shoot(fYp, fSigmaYp); //yp=y'=py/pz=dy/dz
        fGun->SetParticleMomentumDirection(G4ThreeVector(xp, yp, 1.)); //automatically norm.
                    
        //generate the vertex
        fGun->GeneratePrimaryVertex(anEvent);                 
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

