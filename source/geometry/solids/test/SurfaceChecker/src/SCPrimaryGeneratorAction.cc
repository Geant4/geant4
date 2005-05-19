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
// $Id: SCPrimaryGeneratorAction.cc,v 1.1 2005-05-19 13:07:29 link Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Randomize.hh"

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

   G4ThreeVector aSurfacePoint (GetTorusPoint(u,v)) ;  // point A on surface
   G4ThreeVector aNormal(1,0,0) ;  // normal vector

  G4ThreeVector m = (aSurfacePoint-VertexPosition).unit() ;  // direction

  G4double distance = ( aSurfacePoint - VertexPosition ).mag() ;
  G4double theta = std::acos(m.unit()*aNormal) ;

  spoint.SetSurfacePoint(aSurfacePoint) ;
  
  G4PrimaryParticle* aPrimaryParticle =
    new G4PrimaryParticle(aParticleDefinition, m.x(), m.y(), m.z());
  aPrimaryParticle->SetMass (0.);

  aVertex->SetPrimary( aPrimaryParticle );

  G4cout  << "Event "     << i << G4endl 
	  << "Vertex "  << VertexPosition << " " << u << " " << v <<  G4endl 
	  << "Surface "  << aSurfacePoint << G4endl 
	  << "Distance " << distance << G4endl 
	  << "Momentum "  << m << G4endl 
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
  SinTheta = sqrt (1.-CosTheta*CosTheta);
  rand = G4UniformRand();
  Phi = twopi*rand;
  SinPhi = sin (Phi);
  CosPhi = cos (Phi);
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

G4ThreeVector SCPrimaryGeneratorAction::GetTorusPoint(G4double &u, G4double &v)
{

  // get parameters

  G4double c = myDetector->GetTrackerpDx2() ;
  G4double a = myDetector->GetTrackerpDy2() ;

  u = twopi*G4UniformRand();
  v = twopi*G4UniformRand();

  G4ThreeVector retval  ;

  retval.setX((c + a*std::cos(v))* std::cos(u));
  retval.setY((c + a*std::cos(v))* std::sin(u));
  retval.setZ( a * std::sin(v) );

  return retval ;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

