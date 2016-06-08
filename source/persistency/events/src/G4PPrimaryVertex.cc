// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PPrimaryVertex.cc,v 1.7 1999/12/15 14:51:21 gunter Exp $
// GEANT4 tag $Name: geant4-01-01 $
//

#include "G4PPrimaryVertex.hh"
#include <assert.h>
#include "G4ios.hh"

#include "G4PrimaryVertex.hh"
#include "G4PPrimaryParticle.hh"

G4PPrimaryVertex::G4PPrimaryVertex()
:X0(0.),Y0(0.),Z0(0.),T0(0.),numberOfParticle(0),nextVertex(NULL)
// ,theParticle(NULL)
{;}

G4PPrimaryVertex::G4PPrimaryVertex(const G4PrimaryVertex* vertex)
{
  X0 = vertex->GetX0();
  Y0 = vertex->GetY0();
  Z0 = vertex->GetZ0();
  T0 = vertex->GetT0();
  numberOfParticle = vertex->GetNumberOfParticle();

  theParticle = 
         new(ooThis()) G4PPrimaryParticle( vertex->GetPrimary(0) );

  HepRef(G4PPrimaryParticle) particle = theParticle;
  for (G4int i=0; i<numberOfParticle; i++)
  {
    assert( particle != NULL );
    particle = particle->GetNext();
  }
  theTail = particle;

  const G4PrimaryVertex* nextVX = vertex->GetNext();
  if(nextVX)
  { nextVertex = new(ooThis()) G4PPrimaryVertex(nextVX); }
  else
  { nextVertex = NULL; }
}

G4PPrimaryVertex::~G4PPrimaryVertex()
{
  if(theParticle != NULL)
  {
    G4PPrimaryParticle* pp = (HepRef(G4PPrimaryParticle)) theParticle;
    HepDelete(pp);
  }
  { HepDelete(theParticle); }

  if(theTail != NULL)
  {
    G4PPrimaryParticle* pt = (HepRef(G4PPrimaryParticle)) theTail;
    HepDelete(pt);
  }
  { HepDelete(theParticle); }

  if(nextVertex != NULL)
  {
    G4PPrimaryVertex* nv = (HepRef(G4PPrimaryVertex)) nextVertex;
    HepDelete(nv);
  }
}

G4PrimaryVertex* G4PPrimaryVertex::MakeTransientObject()
{
  G4PrimaryVertex* aPV = new G4PrimaryVertex(X0,Y0,Z0,T0);

  if(theParticle != NULL)
  {
    G4PrimaryParticle* particle = theParticle->MakeTransientObject();

    for( G4int i = 0; i<numberOfParticle; i++ )
    {
      assert( particle );
      aPV->SetPrimary( particle );
      particle = particle->GetNext();
    }
  }

  if( nextVertex != NULL )
    aPV->SetNext( nextVertex->MakeTransientObject() );

  return aPV;
}

HepRef(G4PPrimaryParticle) G4PPrimaryVertex::GetPrimary(G4int i=0) const
{
  if( i == 0 )
  { return theParticle; }
  else if( i > 0 && i < numberOfParticle )
  {
    HepRef(G4PPrimaryParticle) particle = theParticle;
    for( int j=0; j<i; j++ )
    {
      if( particle == NULL ) return NULL;
      particle = particle->GetNext();
    }
    return particle;
  }
  else
  { return NULL; }
}

