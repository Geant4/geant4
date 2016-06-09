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
// $Id: FCALPrimaryGeneratorAction.cc,v 1.9 2006-07-21 08:19:33 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FCALPrimaryGeneratorAction.hh"


#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

#include <fstream>
#include <cstdlib>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALPrimaryGeneratorAction::FCALPrimaryGeneratorAction()
{
  X = new G4double[5001];
  Y = new G4double[5001];
  Z = new G4double[5001];
  Cos_X = new G4double[5001];
  Cos_Y = new G4double[5001];
  Cos_Z = new G4double[5001];


  G4int Nparticles = 1;
  particleGun  = new G4ParticleGun(Nparticles);

  // default Particle
  G4String particleName;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();  
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e+");
  particleGun->SetParticleDefinition(particle);
 

  // default Energy
  particleGun->SetParticleEnergy(80*GeV);

  // Read Kinematics from file
  G4int InEvent = 0;
  G4String file_name = "data-tracks/tracks-80GeV.dat";
  std::ifstream Traks_file(file_name);
  if(!Traks_file) G4cerr << "WARNING:  Failed to open file " << file_name << G4endl;
  Traks_file.seekg(0);
  
  while(!(Traks_file.eof())) {
    InEvent++;
    Traks_file >> Ievent >> X[InEvent] >> Y[InEvent] >> Z[InEvent] 
	 >> Cos_X[InEvent] >> Cos_Y[InEvent] >> Cos_Z[InEvent];
  };   

  Nevent = 2500;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALPrimaryGeneratorAction::~FCALPrimaryGeneratorAction()
{
  delete particleGun;
  delete [] X; 
  delete [] Y;
  delete [] Z;
  delete [] Cos_X;
  delete [] Cos_Y;
  delete [] Cos_Z;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
       //this function is called at the begining of event

  Nevent++;
  
  particleGun->SetParticlePosition(G4ThreeVector(X[Nevent]*cm,Y[Nevent]*cm,Z[Nevent]*cm));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(-Cos_X[Nevent],Cos_Y[Nevent],-Cos_Z[Nevent]));
  

  particleGun->GeneratePrimaryVertex(anEvent);

  G4cout << "--------------------------------------------" << G4endl;
  G4cout << " Event,  X,Y,Z Generated Vertex : " << G4endl;
  G4cout << Nevent << " " << X[Nevent] << " " << Y[Nevent] << " " << Z[Nevent]<< G4endl;
  G4cout << "--------------------------------------------" << G4endl;
 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


