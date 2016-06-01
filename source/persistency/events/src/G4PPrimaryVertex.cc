// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PPrimaryVertex.cc,v 2.2 1998/07/13 17:19:19 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//

#include "G4PPrimaryVertex.hh"
#include "G4ios.hh"

G4PPrimaryVertex::G4PPrimaryVertex()
:X0(0.),Y0(0.),Z0(0.),T0(0.),numberOfParticle(0),nextVertex(NULL),
 theParticle(NULL)
{;}

G4PPrimaryVertex::G4PPrimaryVertex(const G4PrimaryVertex* vertex)
{
  X0 = vertex->GetX0();
  Y0 = vertex->GetY0();
  Z0 = vertex->GetZ0();
  T0 = vertex->GetT0();
  numberOfParticle = vertex->GetNumberOfParticle();

  const G4PrimaryVertex* nextVX = vertex->GetNext();
  if(nextVX)
  { nextVertex = new G4PPrimaryVertex(nextVX); }
  else
  { nextVertex = NULL; }

  theParticle = NULL;
}

G4PPrimaryVertex::~G4PPrimaryVertex()
{
  if(theParticle != NULL)
  { delete theParticle; }
  if(nextVertex != NULL)
  { delete nextVertex; }
}



