// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PrimaryVertex.cc,v 1.3 2000-10-19 15:19:37 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4PrimaryVertex.hh"
#include "G4ios.hh"

G4Allocator<G4PrimaryVertex> aPrimaryVertexAllocator;

G4PrimaryVertex::G4PrimaryVertex()
:X0(0.),Y0(0.),Z0(0.),T0(0.),numberOfParticle(0),nextVertex(NULL),
 theParticle(NULL),theTail(NULL),Weight0(1.0)
{;}

G4PrimaryVertex::G4PrimaryVertex(
          G4double x0,G4double y0,G4double z0,G4double t0)
:X0(x0),Y0(y0),Z0(z0),T0(t0),numberOfParticle(0),nextVertex(NULL),
 theParticle(NULL),theTail(NULL),Weight0(1.0)
{;}

G4PrimaryVertex::G4PrimaryVertex(G4ThreeVector xyz0,G4double t0)
:T0(t0),numberOfParticle(0),nextVertex(NULL),
 theParticle(NULL),theTail(NULL),Weight0(1.0)
{
  X0=xyz0.x();
  Y0=xyz0.y();
  Z0=xyz0.z();
}

G4PrimaryVertex::~G4PrimaryVertex()
{
  if(theParticle != NULL)
  { delete theParticle; }
  if(nextVertex != NULL)
  { delete nextVertex; }
}

const G4PrimaryVertex & 
G4PrimaryVertex::operator=(const G4PrimaryVertex &right)
{ return *this; }
int G4PrimaryVertex::operator==(const G4PrimaryVertex &right) const
{ return false; }
int G4PrimaryVertex::operator!=(const G4PrimaryVertex &right) const
{ return true; }

void G4PrimaryVertex::Print() const
{
  G4cout << "Vertex  ( "
       << X0 << ", " << Y0 << ", " << Z0 << ", " << T0 << " )" 
       << " Weight " << Weight0 << G4endl;
  G4cout << "#### Primary particles" << G4endl;
  G4PrimaryParticle* aPrim = theParticle;
  if(aPrim != NULL)
  {
    aPrim->Print();
    aPrim = aPrim->GetNext();
  }
}


