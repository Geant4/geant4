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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//
#include "G4HumanPhantomPrimaryGeneratorAction.hh"
#include "G4HumanPhantomConstruction.hh"
#include "G4HumanPhantomPrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

#include "G4RunManager.hh"
//#ifdef G$ANALYSIS_USE
//#include "G4HumanPhantomAnalysisManager.hh"
//#endif

#include "G4ios.hh"

G4HumanPhantomPrimaryGeneratorAction::G4HumanPhantomPrimaryGeneratorAction()
  :beamKind("beamAlongZ"),worldLength(200.)
{
  G4int n_particle = 1;

  messenger= new G4HumanPhantomPrimaryGeneratorMessenger(this);

  particleGun = new G4ParticleGun(n_particle);

  probability.push_back(0.1667);
  probability.push_back(0.1667);
  probability.push_back(0.1667);
  probability.push_back(0.1667);
  probability.push_back(0.1666);
  probability.push_back(0.1666);

  // Default Particle
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  
particleGun->SetParticleDefinition(particleTable->FindParticle("geantino"));
  particleGun->SetParticleEnergy(100.*MeV);
}

G4HumanPhantomPrimaryGeneratorAction::~G4HumanPhantomPrimaryGeneratorAction()
{
  delete particleGun;
  delete messenger;
}

void G4HumanPhantomPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if (beamKind == "beamAlongZ") GenerateBeamAlongZ();
  if (beamKind == "beamAlongY") GenerateBeamAlongY();
  if (beamKind == "beamAlongX") GenerateBeamAlongX();
  if (beamKind == "isotropicFlux") GenerateIsotropicFlux();

  particleGun->SetParticlePosition(G4ThreeVector(x0, y0,z0));
  particleGun->GeneratePrimaryVertex(anEvent);
}

void G4HumanPhantomPrimaryGeneratorAction::GenerateBeamAlongZ()
{
  z0 = 0.5*(worldLength)*cm;
  y0 = (worldLength)*(G4UniformRand()-0.5)*cm;
  x0 = (worldLength)*(G4UniformRand()-0.5)*cm;
 
  //z0 = 0.5*(worldLength)*cm;
  //y0 = 50. * cm + 5.*(G4UniformRand()-0.5)*cm;
  //x0 = 10.*(G4UniformRand()-0.5)*cm;

   G4ThreeVector direction(0.,0.,-1.);


  particleGun->SetParticleMomentumDirection(direction);
}

void G4HumanPhantomPrimaryGeneratorAction::GenerateBeamAlongX()
{
  x0 = -0.5*(worldLength)*cm;
  y0 = (worldLength)*(G4UniformRand()-0.5)*cm;
  z0 = (worldLength)*(G4UniformRand()-0.5)*cm;
 
  G4ThreeVector direction(1.,0.,0.);
  particleGun->SetParticleMomentumDirection(direction);
}
void G4HumanPhantomPrimaryGeneratorAction::GenerateBeamAlongY()
{
  y0 = 0.5*(worldLength)*cm;
  x0 = (worldLength)*(G4UniformRand()-0.5)*cm;
  z0 = (worldLength)*(G4UniformRand()-0.5)*cm;
 
  G4ThreeVector direction(0.,-1.,0.);
  particleGun->SetParticleMomentumDirection(direction);
}

void G4HumanPhantomPrimaryGeneratorAction::GenerateIsotropicFlux()
{
  G4double random = G4UniformRand();
  G4double sum = 0.;
  G4int i = 0;

  while(sum<random){sum += probability[i]; i++;}
  
  if(i==1) 
    {
      z0 = -0.5*(worldLength-2.)*cm;
      y0 = (worldLength)*(G4UniformRand()-0.5)*cm;
      x0 = (worldLength)*(G4UniformRand()-0.5)*cm;
    }

  if(i==2) 
    {
      y0 = -0.5*(worldLength-2.)*cm;
      z0 = (worldLength)*(G4UniformRand()-0.5)*cm;
      x0 = (worldLength)*(G4UniformRand()-0.5)*cm;
    }

  if(i==3) 
    {
      x0 = -0.5*(worldLength-2.)*cm;
      z0 = (worldLength)*(G4UniformRand()-0.5)*cm;
      y0 = (worldLength)*(G4UniformRand()-0.5)*cm;
    }

  if (i==4)
    {
      z0 = 0.5*(worldLength-2.)*cm;
      y0 = (worldLength)*(G4UniformRand()-0.5)*cm;
      x0 = (worldLength)*(G4UniformRand()-0.5)*cm;
    }

 if(i==5) 
    {
      y0 = 0.5*(worldLength-2.)*cm;
      z0 = (worldLength)*(G4UniformRand()-0.5)*cm;
      x0 = (worldLength)*(G4UniformRand()-0.5)*cm;
    }

  if(i==6) 
    {
      x0 = 0.5*(worldLength-2.)*cm;
      z0 = (worldLength)*(G4UniformRand()-0.5)*cm;
      y0 = (worldLength)*(G4UniformRand()-0.5)*cm;
    } 

  G4double a,b,c;
  G4double n;
  do{
    a = (G4UniformRand()-0.5)/0.5;
    b = (G4UniformRand()-0.5)/0.5; 
    c = (G4UniformRand()-0.5)/0.5;
    n = a*a+b*b+c*c;
  }while(n > 1 || n == 0.0);
  n = std::sqrt(n);
  a /= n;
  b /= n;
  c /= n;

  G4ThreeVector direction(a,b,c);
  particleGun->SetParticleMomentumDirection(direction);  
}

void G4HumanPhantomPrimaryGeneratorAction::SetBeam(G4String beam)
{
  if((beam == "beamAlongZ")||(beam == "beamAlongX")||
  (beam == "beamAlongY")||(beam == "isotropicFlux")) beamKind = beam;
  
  else G4cout<<"This option is not valid "<<
               "---> beamAlongZ/beamAlongY/beamAlongX/isotropicFlux"
             <<G4endl;
}

