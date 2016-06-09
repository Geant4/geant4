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
// Code developed by:
//  S.Larsson
//
//    ********************************************
//    *                                          *
//    *    PurgMagPrimaryGeneratorAction.cc     *
//    *                                          *
//    ********************************************
//
// $Id: PurgMagPrimaryGeneratorAction.cc,v 1.2 2004/06/18 09:17:59 gunter Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "PurgMagPrimaryGeneratorAction.hh"

#include "PurgMagDetectorConstruction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"

//Print position of primaries.
#define POSITION 0

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PurgMagPrimaryGeneratorAction::PurgMagPrimaryGeneratorAction(PurgMagDetectorConstruction* PurgMagDC)
  :PurgMagDetector(PurgMagDC),rndmVertex(false)
{
  //default kinematic
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  G4ParticleDefinition* particle
    =
    
    G4ParticleTable::GetParticleTable()->FindParticle("e-");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleEnergy(50.*MeV);

  //Momentum Direction
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PurgMagPrimaryGeneratorAction::~PurgMagPrimaryGeneratorAction()
{
  delete particleGun;
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PurgMagPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  //Start position of primaries
  G4double z0 = 15.*cm;
  G4double x0 = 0.*cm;
  G4double y0 = 0.*cm;
  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  particleGun->GeneratePrimaryVertex(anEvent);

#if POSITION  
  G4cout << "\n----Particle Gun--------------------------------------------\n";
  G4cout << "\n ---> SetParticlePosition(G4ThreeVector(x0,y0,z0)) ";
  G4cout << "\n ---> x0 = " << x0 << " "<< G4BestUnit(x0,"Length");
  G4cout << "\n ---> y0 = " << y0 << " "<< G4BestUnit(y0,"Length");
  G4cout << "\n ---> z0 = " << z0 << " "<< G4BestUnit(z0,"Length");
  G4cout << "\n-----------------------------------------------------------\n";
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


