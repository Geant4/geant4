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
/// \file electromagnetic/TestEm10/src/Em10PrimaryGeneratorAction.cc
/// \brief Implementation of the Em10PrimaryGeneratorAction class
//
//
// $Id: Em10PrimaryGeneratorAction.cc 73033 2013-08-15 09:24:45Z gcosmo $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em10PrimaryGeneratorAction.hh"

#include "Em10DetectorConstruction.hh"
#include "Em10PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

 G4String Em10PrimaryGeneratorAction::thePrimaryParticleName="proton";

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em10PrimaryGeneratorAction::Em10PrimaryGeneratorAction(
                                            Em10DetectorConstruction*)
//                                            Em10DetectorConstruction* Em10DC)
:G4VUserPrimaryGeneratorAction(),
// Em10Detector(Em10DC),
 rndmFlag("off"),xvertex(0.),yvertex(0.),zvertex(0.),
 vertexdefined(false)
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);

  //create a messenger for this class
  gunMessenger = new Em10PrimaryGeneratorMessenger(this);

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
                    = particleTable->FindParticle(particleName="proton");
  particleGun->SetParticleDefinition(particle);

  thePrimaryParticleName = particle->GetParticleName();

  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticleEnergy(100.*GeV);

  zvertex = 0.0 ; //  -0.5*(Em10Detector->GetAbsorberThickness());
  particleGun->SetParticlePosition(G4ThreeVector(xvertex,yvertex,zvertex));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em10PrimaryGeneratorAction::~Em10PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  //
  thePrimaryParticleName = particleGun->GetParticleDefinition()->
                                                GetParticleName();
  /* ****************************************************
  G4double x0,y0,z0;
  if(vertexdefined)
  {
    x0 = xvertex;
    y0 = yvertex;
    z0 = zvertex;
  }
  else
  {
    x0 = 0.;
    y0 = 0.;
    z0 = 0.; // -0.5*(Em10Detector->GetWorldSizeZ()) ;
  }
  G4double r0,phi0;
  if (rndmFlag == "on")
  {
      r0 = Em10Detector->GetAbsorberRadius())*std::sqrt(G4UniformRand();
      phi0 = twopi*G4UniformRand();
      x0 = r0*std::cos(phi0);
      y0 = r0*std::sin(phi0);
  }
  ********************************************* */
  //  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  particleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
  particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String Em10PrimaryGeneratorAction::GetPrimaryName()
{
   return thePrimaryParticleName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10PrimaryGeneratorAction::Setzvertex(G4double z)
{
  vertexdefined = true;
  zvertex = z;
  G4cout << " Z coordinate of the primary vertex = " << zvertex/mm <<
            " mm." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10PrimaryGeneratorAction::Setxvertex(G4double x)
{
  vertexdefined = true;
  xvertex = x;
  G4cout << " X coordinate of the primary vertex = " << xvertex/mm <<
            " mm." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10PrimaryGeneratorAction::Setyvertex(G4double y)
{
  vertexdefined = true;
  yvertex = y;
  G4cout << " Y coordinate of the primary vertex = " << yvertex/mm <<
            " mm." << G4endl;
}
