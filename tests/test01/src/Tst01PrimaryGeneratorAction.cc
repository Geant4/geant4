
//
//
//
//

#include "Tst01PrimaryGeneratorAction.hh"

#include "Tst01PrimaryGeneratorMessenger.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4TransportationManager.hh"
#include "globals.hh"
#include "CLHEP/Random/RandFlat.h"

////////////////////////////////////////////////////////////////////////
//
//

Tst01PrimaryGeneratorAction::Tst01PrimaryGeneratorAction():
  generatorAction (standardGun),
  fGunPosition(0.0, 0.0, 0.0),
  particleGun (0),
  messenger (0),
  worldVolume (0)
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  // default particle

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle 
    = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  particleGun->SetParticleEnergy(1.*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));

  // messenger

  messenger = new Tst01PrimaryGeneratorMessenger (this);

  // world extent

  worldVolume = G4TransportationManager::GetTransportationManager ()
              -> GetNavigatorForTracking () -> GetWorldVolume ()     ;
  if (worldVolume) 
  {  
    worldExtent = worldVolume -> GetLogicalVolume () -> 
                                         GetSolid () -> GetExtent ();

    fSize = sqrt( ( worldExtent.GetXmax() - worldExtent.GetXmin() )*
                  ( worldExtent.GetXmax() - worldExtent.GetXmin() ) +

                  ( worldExtent.GetYmax() - worldExtent.GetYmin() )*
                  ( worldExtent.GetYmax() - worldExtent.GetYmin() ) +

                  ( worldExtent.GetZmax() - worldExtent.GetZmin() )* 
                  ( worldExtent.GetZmax() - worldExtent.GetZmin() )     ) ;
  }
  else
  {
    fSize = 0.0 ;
  }
  G4cout<<"fSize = "<<fSize<<endl ;
}

///////////////////////////////////////////////////////////////////////////
//
// Destructor: delets pointers

Tst01PrimaryGeneratorAction::~Tst01PrimaryGeneratorAction()
{
  delete particleGun;
  delete messenger;
}

////////////////////////////////////////////////////////////////////////
//
//

void Tst01PrimaryGeneratorAction::SelectPrimaryGeneratorAction
(Action action)
{
  generatorAction = action;
}

////////////////////////////////////////////////////////////////////////
//
//

void Tst01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4VPhysicalVolume* currentWorldVolume;
  G4ThreeVector direction, position ;
  G4double costheta, sintheta, phi, cosphi, sinphi ;
  G4double rho;

  switch (generatorAction) 
  {

    case standardGun:  

    particleGun->GeneratePrimaryVertex(anEvent);
    break;
//...................................................................
    case randomDirectionGun:   

    costheta = RandFlat::shoot (-1., 1.);
    sintheta = sqrt (1. - costheta * costheta);
    phi      = RandFlat::shoot (twopi);
    cosphi   = cos (phi);
    sinphi   = sin (phi);
    particleGun->SetParticleMomentumDirection
      (G4ThreeVector (sintheta * cosphi, sintheta * sinphi, costheta));

    particleGun->GeneratePrimaryVertex(anEvent);
    break;
//......................................................................
    case randomPositionGun:   

    // Check if world is in place or has changed.
    currentWorldVolume = 
      G4TransportationManager::GetTransportationManager ()
      -> GetNavigatorForTracking () -> GetWorldVolume ();

    if (!worldVolume ||	worldVolume != currentWorldVolume) 
    {
      worldVolume = currentWorldVolume;
      if (worldVolume) worldExtent = worldVolume -> GetLogicalVolume ()
			 -> GetSolid () -> GetExtent ();
    }

    particleGun->SetParticlePosition
      (G4ThreeVector
       (RandFlat::shoot (worldExtent.GetXmin (), worldExtent.GetXmax ()),
	RandFlat::shoot (worldExtent.GetYmin (), worldExtent.GetYmax ()),
	RandFlat::shoot (worldExtent.GetZmin (), worldExtent.GetZmax ())));

    particleGun->GeneratePrimaryVertex(anEvent);
    break;
//.................................................................
    case randomPositionAndDirectionGun:  

    // Check if world is in place or has changed.
    currentWorldVolume =
      G4TransportationManager::GetTransportationManager ()
      -> GetNavigatorForTracking () -> GetWorldVolume ();

    if (!worldVolume ||	worldVolume != currentWorldVolume) 
    {
      worldVolume = currentWorldVolume;
      if (worldVolume) worldExtent = worldVolume -> GetLogicalVolume ()
			 -> GetSolid () -> GetExtent ();
    }

    particleGun->SetParticlePosition
      (G4ThreeVector
       (RandFlat::shoot (worldExtent.GetXmin (), worldExtent.GetXmax ()),
	RandFlat::shoot (worldExtent.GetYmin (), worldExtent.GetYmax ()),
	RandFlat::shoot (worldExtent.GetZmin (), worldExtent.GetZmax ())));

    costheta = RandFlat::shoot (-1., 1.);
    sintheta = sqrt (1. - costheta * costheta);
    phi      = RandFlat::shoot (twopi);
    cosphi   = cos (phi);
    sinphi   = sin (phi);
    particleGun->SetParticleMomentumDirection
      (G4ThreeVector (sintheta * cosphi, sintheta * sinphi, costheta));

    particleGun->GeneratePrimaryVertex(anEvent);
    break;
//.................................................................
    case viewerGun:

    particleGun->SetParticlePosition(fGunPosition) ;
    if(fPosition)
    {
      direction = G4ThreeVector( -fGunPosition.x()/fPosition ,
                                 -fGunPosition.y()/fPosition ,
                                 -fGunPosition.z()/fPosition    ) ;

      costheta = RandFlat::shoot (fPosition/sqrt(fPosition*fPosition +
                                                  fSize*fSize), 1.);
      sintheta = sqrt (1. - costheta * costheta);
      phi      = RandFlat::shoot (twopi);
      cosphi   = cos (phi) ;
      sinphi   = sin (phi) ;
      position = G4ThreeVector (sintheta * cosphi, 
                                sintheta * sinphi, costheta) ;

      position.rotateUz(direction) ;

      particleGun->SetParticleMomentumDirection(position);      
    }
    else
    {
      costheta = RandFlat::shoot (-1., 1.);
      sintheta = sqrt (1. - costheta * costheta);
      phi      = RandFlat::shoot (twopi);
      cosphi   = cos (phi);
      sinphi   = sin (phi);
      particleGun->SetParticleMomentumDirection
      (G4ThreeVector (sintheta * cosphi, sintheta * sinphi, costheta));
    }
    particleGun->GeneratePrimaryVertex(anEvent) ;
    break ;
//..................................................................
    case planeGun:

    if(fPosition)
    {
      direction = G4ThreeVector( -fGunPosition.x()/fPosition ,
                                 -fGunPosition.y()/fPosition ,
                                 -fGunPosition.z()/fPosition    ) ;
      particleGun->SetParticleMomentumDirection(direction) ;

   // rho = fSize*randFlat::shoot(0.0,1.0) ; phi = RandFlat::shoot (twopi) ; 
   // cosphi   = cos (phi); sinphi   = sin (phi);
   // position = G4ThreeVector(rho*cosphi,rho*sinphi,0.0) ;

      position = G4ThreeVector(RandFlat::shoot (-0.5*fSize,0.5*fSize) ,
                               RandFlat::shoot (-0.5*fSize,0.5*fSize) , 0.0 ) ;
      position.rotateUz(direction) ;
      particleGun->SetParticlePosition(fGunPosition+position) ;
    }
    else
    {
      G4Exception("Invalid setting for plane gun in Tst01PrimaryGeneratorAction::GeneratePrimaries") ;
    }
    particleGun->GeneratePrimaryVertex(anEvent);
    break ;
  }
}


//
//
/////////////////////////////////////////////////////////////////////////
