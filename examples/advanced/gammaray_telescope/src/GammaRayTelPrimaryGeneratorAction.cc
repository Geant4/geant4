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
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelPrimaryGeneratorAction  ------
//           by  G.Santin, F.Longo & R.Giannitrapani (13 nov 2000)
//
// ************************************************************

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4RunManager.hh"
#include "GammaRayTelPrimaryGeneratorAction.hh"

#include "GammaRayTelDetectorConstruction.hh"
#include "GammaRayTelPrimaryGeneratorMessenger.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelPrimaryGeneratorAction::GammaRayTelPrimaryGeneratorAction()
//  :rndmFlag("off"),nSourceType(0),nSpectrumType(0),sourceGun(false)
{
  GammaRayTelDetector = static_cast<const GammaRayTelDetectorConstruction*>
    (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  //create a messenger for this class
  
  gunMessenger = new GammaRayTelPrimaryGeneratorMessenger(this);

  rndmFlag = "off";
  nSourceType = 0;
  nSpectrumType = 0;
  sourceGun = false;
  dVertexRadius = 15.*cm; 
  
  G4int n_particle = 1;

  particleGun  = new G4ParticleGun(n_particle);     
  // default particle kinematic
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  particleGun->SetParticleEnergy(30.*MeV);
  G4double position = 0.5*(GammaRayTelDetector->GetWorldSizeZ());
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,position));
  particleSource = new G4GeneralParticleSource();
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelPrimaryGeneratorAction::~GammaRayTelPrimaryGeneratorAction()
{
 
  delete particleGun;
  delete particleSource;
  delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
   if (sourceGun)
    {

      G4cout << "Using G4ParticleGun ... " << G4endl; 

      //this function is called at the begining of event
      // 
      G4double z0 = 0.5*(GammaRayTelDetector->GetWorldSizeZ());
      G4double x0 = 0.*cm, y0 = 0.*cm;
  
      G4ThreeVector pos0;
      G4ThreeVector vertex0 = G4ThreeVector(x0,y0,z0);
      G4ThreeVector dir0 = G4ThreeVector(0.,0.,-1.);
      
      G4double theta, phi;
      G4double y = 0.;
      G4double f = 0.;
      G4double theta0=0.;
      G4double phi0=0.;
      
      switch(nSourceType) {
      case 0:
	particleGun->SetParticlePosition(vertex0);
	particleGun->SetParticleMomentumDirection(dir0);
	break;
      case 1:
	// GS: Generate random position on the 4PIsphere to create a unif. distrib.
	// GS: on the sphere
	phi = G4UniformRand() * twopi;
	do {
	  y = G4UniformRand()*1.0;
	  theta = G4UniformRand() * pi;
	  f = std::sin(theta);
	} while (y > f);
	vertex0 = G4ThreeVector(1.,0.,0.);
	vertex0.setMag(dVertexRadius);
	vertex0.setTheta(theta);
	vertex0.setPhi(phi);
	particleGun->SetParticlePosition(vertex0);
	
	dir0 = G4ThreeVector(1.,0.,0.);
	do {
	  phi = G4UniformRand() * twopi;
	  do {
	    y = G4UniformRand()*1.0;
	    theta = G4UniformRand() * pi;
	    f = std::sin(theta);
	  } while (y > f);
	  dir0.setPhi(phi);
	  dir0.setTheta(theta);
	} while (vertex0.dot(dir0) >= -0.7 * vertex0.mag());
	particleGun->SetParticleMomentumDirection((G4ParticleMomentum)dir0);
	
	break;
      case 2:
	// GS: Generate random position on the upper semi-sphere z>0 to create a unif. distrib.
	// GS: on a plane
	phi = G4UniformRand() * twopi;
	do {
	  y = G4UniformRand()*1.0;
	  theta = G4UniformRand() * halfpi;
	  f = std::sin(theta) * std::cos(theta);
	} while (y > f);
	vertex0 = G4ThreeVector(1.,0.,0.);
	
	G4double xy = GammaRayTelDetector->GetWorldSizeXY();
	G4double z = GammaRayTelDetector->GetWorldSizeZ();
	
	if (dVertexRadius > xy*0.5)
	  { 
	    G4cout << "vertexRadius too big " << G4endl;
	    G4cout << "vertexRadius setted to " << xy*0.45 << G4endl;
	    dVertexRadius = xy*0.45;
	  }
	
	if (dVertexRadius > z*0.5)
	  { 
	    G4cout << "vertexRadius too high " << G4endl;
	    G4cout << "vertexRadius setted to " << z*0.45 << G4endl;
	    dVertexRadius = z*0.45;
	  }
	
	
	vertex0.setMag(dVertexRadius);
	vertex0.setTheta(theta);
	vertex0.setPhi(phi);
	
	// GS: Get the user defined direction for the primaries and
	// GS: Rotate the random position according to the user defined direction for the particle
	
	dir0 = particleGun->GetParticleMomentumDirection();
	if (dir0.mag() > 0.001) 
	  {
	    theta0 = dir0.theta();
	    phi0   = dir0.phi();   
	  }
	
	if (theta0!=0.) 
	  {
	    G4ThreeVector rotationAxis(1.,0.,0.);
	    rotationAxis.setPhi(phi0+halfpi);
	    vertex0.rotate(theta0+pi,rotationAxis);
	  }
	particleGun->SetParticlePosition(vertex0);
	break;
      }
      
      
      G4double pEnergy = 100*MeV;
      
      switch(nSpectrumType) {
      case 0: // Uniform energy (1-10 GeV)
	y = G4UniformRand();
	pEnergy = y*9.0*GeV + 1.0*GeV;
	G4cout << pEnergy/GeV << " LIN" << G4endl;
	break;
      case 1: // Logaritmic energy
	y = G4UniformRand();
	pEnergy = std::pow(10,y)*GeV;
	G4cout << pEnergy/GeV << " LOG" << G4endl;
	break;
      case 2: // Power Law (-4)
	do {
	  y = G4UniformRand()*100000.0;
	  pEnergy = G4UniformRand() * 10. * GeV;
	  f = std::pow(pEnergy * (1/GeV), -4.);
	} while (y > f);
	//	particleGun->SetParticleEnergy(pEnergy);
	break;
      case 3: // Monochromatic 
	pEnergy = particleGun->GetParticleEnergy(); 
	//100 * MeV; 
	G4cout << pEnergy << " MONO" << G4endl;
	break;
      }
      particleGun->SetParticleEnergy(pEnergy);
      G4cout << particleGun->GetParticleDefinition()->GetParticleName() << G4endl;
      particleGun->GeneratePrimaryVertex(anEvent);
    }
   else
     {
       particleSource->GeneratePrimaryVertex(anEvent);
     }
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....








