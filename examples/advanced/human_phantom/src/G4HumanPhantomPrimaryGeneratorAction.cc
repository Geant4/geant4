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
  particleGun->SetParticleDefinition(particleTable->FindParticle("proton"));
  particleGun->SetParticleEnergy(100.*MeV);
}

G4HumanPhantomPrimaryGeneratorAction::~G4HumanPhantomPrimaryGeneratorAction()
{
  delete particleGun;
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
  n = sqrt(n);
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

