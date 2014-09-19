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
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "SCPrimaryGeneratorAction.hh"
#include "SCDetectorConstruction.hh"

#include "SCSurfacePoint.hh" 
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SurfacePoint spoint ;

SCPrimaryGeneratorAction::SCPrimaryGeneratorAction(
                                               SCDetectorConstruction* myDC)
:myDetector(myDC)
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

// default particle

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("geantino");
  
  particleGun->SetParticleDefinition(particle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SCPrimaryGeneratorAction::~SCPrimaryGeneratorAction()
{
  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* aParticleDefinition 
    = particleTable->FindParticle(particleName="geantino");

  // create a new vertex in a random position 

  G4int i = anEvent->GetEventID() ;

  G4ThreeVector VertexPosition (GetRandomPosition());
  G4PrimaryVertex* aVertex = new G4PrimaryVertex( VertexPosition, 0);

  G4double u,v ;

  G4ThreeVector aSurfacePoint (GetSurfacePoint(u,v)) ;  // point A on surface
  G4ThreeVector aNormal(1,0,0) ;  // normal vector

  G4ThreeVector md = (aSurfacePoint-VertexPosition).unit() ;  // direction

  G4double distance = ( aSurfacePoint - VertexPosition ).mag() ;
  G4double theta = std::acos(md.unit()*aNormal) ;

  spoint.SetSurfacePoint(aSurfacePoint) ;
  
  G4PrimaryParticle* aPrimaryParticle =
    new G4PrimaryParticle(aParticleDefinition, md.x(), md.y(), md.z());
  aPrimaryParticle->SetMass (0.);

  aVertex->SetPrimary( aPrimaryParticle );

  G4cout  << "Event "     << i << G4endl 
	  << "Vertex "  << VertexPosition << " " << u << " " << v <<  G4endl 
	  << "Surface "  << aSurfacePoint << G4endl 
	  << "Distance " << distance << G4endl 
	  << "Momentum "  << md << G4endl 
          << "Angle " << theta 
	  << G4endl ;

  anEvent->AddPrimaryVertex( aVertex );


}


G4ThreeVector SCPrimaryGeneratorAction::GetRandomDirection() {

  G4ThreeVector retval;

  G4double CosTheta;
  G4double SinTheta;

  G4double Phi;
  G4double SinPhi;
  G4double CosPhi;

  G4double rand;

  rand = G4UniformRand();

  CosTheta = 2.0*rand -1.0;
  SinTheta = std::sqrt (1.-CosTheta*CosTheta);
  rand = G4UniformRand();
  Phi = twopi*rand;
  SinPhi = std::sin (Phi);
  CosPhi = std::cos (Phi);
  retval.setX(SinTheta*CosPhi);
  retval.setY(SinTheta*SinPhi);
  retval.setZ(CosTheta);

  return retval;
}

G4ThreeVector SCPrimaryGeneratorAction::GetRandomPosition() 
{

  G4double a = 0.5*myDetector->GetWorldFullLength();

  G4double x = ( G4UniformRand()*2 - 1 )*a;
  G4double y = ( G4UniformRand()*2 - 1 )*a;
  G4double z = ( G4UniformRand()*2 - 1 )*a;

  G4ThreeVector retval (x, y, z);

  return retval;
}

G4ThreeVector SCPrimaryGeneratorAction::GetSurfacePoint(G4double &, G4double &)
{


  G4String val = myDetector->GetDetectorType() ;


  return myDetector->GetSolid()->GetPointOnSurface() ;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

