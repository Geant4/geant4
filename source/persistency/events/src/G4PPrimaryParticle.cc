// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PPrimaryParticle.cc,v 1.1 1999/12/05 22:32:26 morita Exp $
// GEANT4 tag $Name: geant4-03-01 $
//

#include "G4PPrimaryParticle.hh"

#include "G4PrimaryParticle.hh"

G4PPrimaryParticle::G4PPrimaryParticle(G4PrimaryParticle* aParticle)
{
  PDGcode = aParticle->GetPDGcode();

  G4ParticleDefinition* ag4code  = aParticle->GetG4code();
  // if( ag4code != NULL )
  // { ...will look up G4PParticleDefinition and set the G4code...; }

  Px      = aParticle->GetPx(); 
  Py      = aParticle->GetPy(); 
  Pz      = aParticle->GetPz(); 

  G4PrimaryParticle* aNext = aParticle->GetNext();
  if( aNext != NULL )
  { nextParticle = new(ooThis()) G4PPrimaryParticle( aNext ); }

  G4PrimaryParticle* aDaughter = aParticle->GetDaughter();
  if( aDaughter != NULL )
  { daughterParticle = new(ooThis()) G4PPrimaryParticle( aDaughter ); }

  trackID = aParticle->GetTrackID();
  mass    = aParticle->GetMass();
  polX    = aParticle->GetPolX();
  polY    = aParticle->GetPolY();
  polZ    = aParticle->GetPolZ();
}

G4PPrimaryParticle::~G4PPrimaryParticle()
{
  if(nextParticle != NULL)
  { HepDelete(nextParticle); }
  if(daughterParticle != NULL)
  { HepDelete(daughterParticle); }
}

G4PrimaryParticle* G4PPrimaryParticle::MakeTransientObject()
{
  G4PrimaryParticle* aParticle = new G4PrimaryParticle(PDGcode, Px, Py, Pz);

  // aParticle->SetG4Code(...will lookup G4ParticleDefinition and set it...;)

  if(nextParticle != NULL)
  { aParticle->SetNext( nextParticle->MakeTransientObject() ); }
  if(daughterParticle != NULL)
  { aParticle->SetDaughter( daughterParticle->MakeTransientObject() ); }

  aParticle->SetTrackID( trackID );
  aParticle->SetMass( mass );
  aParticle->SetPolarization( polX, polY, polZ );

  return aParticle;
}

