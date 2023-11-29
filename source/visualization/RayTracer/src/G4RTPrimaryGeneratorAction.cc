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

#include "G4RTPrimaryGeneratorAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4TransportationManager.hh"
#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"

#include "G4TheMTRayTracer.hh"

G4RTPrimaryGeneratorAction::G4RTPrimaryGeneratorAction()
{
  G4ThreeVector zero;
  particle_definition = 0;
  particle_energy = 1.0*CLHEP::GeV;
  particle_time = 0.0;
  particle_polarization = zero;

  pWorld = 0;
  whereisit = kInside;

  nRow = 0;
  nColumn = 0;
  
  eyePosition = zero;
  eyeDirection = zero;
  up = G4ThreeVector(0,1,0);
  headAngle = 0.0;
  viewSpan = 0.0;
  stepAngle = 0.0;
  viewSpanX = 0.0;
  viewSpanY = 0.0;

  distortionOn = false;
}

G4RTPrimaryGeneratorAction::~G4RTPrimaryGeneratorAction()
{;}

void G4RTPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // Note: We don't use G4ParticleGun here, as instantiating a G4ParticleGun
  //  object causes creation of UI commands and corresponding UI messenger
  //  that interfare with normal G4ParticleGun UI commands.

  // evId = iRow * nColumn + iColumn
  G4int evId = anEvent->GetEventID(); 
  G4int iRow = evId / nColumn;
  G4int iColumn = evId % nColumn;
  G4double angleX = -(viewSpanX/2. - G4double(iColumn)*stepAngle);
  G4double angleY = viewSpanY/2. - G4double(iRow)*stepAngle;
  G4ThreeVector rayDirection;
  if(distortionOn)
  { rayDirection = G4ThreeVector(-std::tan(angleX)/std::cos(angleY),std::tan(angleY)/std::cos(angleX),1.0); }
  else
  { rayDirection = G4ThreeVector(-std::tan(angleX),std::tan(angleY),1.0); }
  G4double cp = std::cos(eyeDirection.phi());
  G4double sp = std::sqrt(1.-cp*cp);
  G4double ct = std::cos(eyeDirection.theta());
  G4double st = std::sqrt(1.-ct*ct);
  G4double gam = std::atan2(ct*cp*up.x()+ct*sp*up.y()-st*up.z(), -sp*up.x()+cp*up.y());
  rayDirection.rotateZ(-gam);
  rayDirection.rotateZ(headAngle);
  rayDirection.rotateUz(eyeDirection);

  G4ThreeVector rayPosition(eyePosition);
  if (whereisit != kInside) {
    // Eye position is outside the world, so move it inside.
    G4double outsideDistance = pWorld->GetLogicalVolume()->GetSolid()->
       DistanceToIn(rayPosition,rayDirection);
    if(outsideDistance != kInfinity)
    { rayPosition = rayPosition + (outsideDistance+0.001)*rayDirection; }
    else
    {
      // Ray does not intercept world at all.
      // Return without primary particle.
      return;
    }
  }
  
  // create a new vertex
  G4PrimaryVertex* vertex = new G4PrimaryVertex(rayPosition,particle_time);

  // create new primaries and set them to the vertex
  G4double mass = particle_definition->GetPDGMass();
  G4PrimaryParticle* particle = new G4PrimaryParticle(particle_definition);
  particle->SetKineticEnergy( particle_energy );
  particle->SetMass( mass );
  particle->SetMomentumDirection( rayDirection.unit() );
  particle->SetPolarization(particle_polarization.x(),
                            particle_polarization.y(),
                            particle_polarization.z());
  vertex->SetPrimary( particle );

  anEvent->AddPrimaryVertex( vertex );
}

void G4RTPrimaryGeneratorAction::SetUp()
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  particle_definition = particleTable->FindParticle("geantino");
  if(!particle_definition)
  {
    G4String msg;
    msg =  " G4RayTracer uses geantino to trace the ray, but your physics list does not\n";
    msg += "define G4Geantino. Please add G4Geantino in your physics list.";
    G4Exception("G4RTPrimaryGeneratorAction::SetUp","VisRayTracer00101",FatalException,msg);
  }

  G4TheMTRayTracer* rt = G4TheMTRayTracer::theInstance;
  nRow = rt->nRow;
  nColumn = rt->nColumn;
  eyePosition = rt->eyePosition;
  eyeDirection = rt->eyeDirection;
  viewSpan = rt->viewSpan;
  stepAngle = viewSpan/100.;
  viewSpanX = stepAngle*nColumn;
  viewSpanY = stepAngle*nRow;
  distortionOn = rt->distortionOn;

  pWorld = G4TransportationManager::GetTransportationManager()->
        GetNavigatorForTracking()->GetWorldVolume();
  whereisit = pWorld->GetLogicalVolume()->GetSolid()->Inside(eyePosition);
}

