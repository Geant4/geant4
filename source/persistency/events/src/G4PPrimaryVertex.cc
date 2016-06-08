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
// $Id: G4PPrimaryVertex.cc,v 1.9 2001/07/11 10:02:16 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $
//

#include "G4PPrimaryVertex.hh"
#include <assert.h>
#include "G4ios.hh"

#include "G4PrimaryVertex.hh"
#include "G4PPrimaryParticle.hh"

G4PPrimaryVertex::G4PPrimaryVertex()
:X0(0.),Y0(0.),Z0(0.),T0(0.),numberOfParticle(0),nextVertex(0)
// ,theParticle(0)
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
    assert( particle != 0 );
    particle = particle->GetNext();
  }
  theTail = particle;

  const G4PrimaryVertex* nextVX = vertex->GetNext();
  if(nextVX)
  { nextVertex = new(ooThis()) G4PPrimaryVertex(nextVX); }
  else
  { nextVertex = 0; }
}

G4PPrimaryVertex::~G4PPrimaryVertex()
{
  if(theParticle != 0)
  {
    G4PPrimaryParticle* pp = (HepRef(G4PPrimaryParticle)) theParticle;
    HepDelete(pp);
  }
  { HepDelete(theParticle); }

  if(theTail != 0)
  {
    G4PPrimaryParticle* pt = (HepRef(G4PPrimaryParticle)) theTail;
    HepDelete(pt);
  }
  { HepDelete(theParticle); }

  if(nextVertex != 0)
  {
    G4PPrimaryVertex* nv = (HepRef(G4PPrimaryVertex)) nextVertex;
    HepDelete(nv);
  }
}

G4PrimaryVertex* G4PPrimaryVertex::MakeTransientObject()
{
  G4PrimaryVertex* aPV = new G4PrimaryVertex(X0,Y0,Z0,T0);

  if(theParticle != 0)
  {
    G4PrimaryParticle* particle = theParticle->MakeTransientObject();

    for( G4int i = 0; i<numberOfParticle; i++ )
    {
      assert( particle );
      aPV->SetPrimary( particle );
      particle = particle->GetNext();
    }
  }

  if( nextVertex != 0 )
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
      if( particle == 0 ) return 0;
      particle = particle->GetNext();
    }
    return particle;
  }
  else
  { return 0; }
}

