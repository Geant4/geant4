// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PrimaryParticle.cc,v 1.1 1999-01-07 16:06:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4PrimaryParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ios.hh"

G4Allocator<G4PrimaryParticle> aPrimaryParticleAllocator;

G4PrimaryParticle::G4PrimaryParticle()
:PDGcode(0),G4code(NULL),Px(0.),Py(0.),Pz(0.),
 nextParticle(NULL),daughterParticle(NULL),trackID(-1),
 mass(0.),polX(0.),polY(0.),polZ(0.)
{;}

G4PrimaryParticle::G4PrimaryParticle(G4int Pcode)
:PDGcode(Pcode),Px(0.),Py(0.),Pz(0.),
 nextParticle(NULL),daughterParticle(NULL),trackID(-1),
 mass(0.),polX(0.),polY(0.),polZ(0.)
{ G4code = G4ParticleTable::GetParticleTable()->FindParticle(Pcode); }

G4PrimaryParticle::G4PrimaryParticle(G4int Pcode,
                        G4double px,G4double py,G4double pz)
:PDGcode(Pcode),Px(px),Py(py),Pz(pz),
 nextParticle(NULL),daughterParticle(NULL),trackID(-1),
 mass(0.),polX(0.),polY(0.),polZ(0.)
{ G4code = G4ParticleTable::GetParticleTable()->FindParticle(Pcode); }

G4PrimaryParticle::G4PrimaryParticle(G4ParticleDefinition* Gcode)
:G4code(Gcode),Px(0.),Py(0.),Pz(0.),
 nextParticle(NULL),daughterParticle(NULL),trackID(-1),
 mass(0.),polX(0.),polY(0.),polZ(0.)
{ PDGcode = Gcode->GetPDGEncoding(); }

G4PrimaryParticle::G4PrimaryParticle(G4ParticleDefinition* Gcode,
                        G4double px,G4double py,G4double pz)
:G4code(Gcode),Px(px),Py(py),Pz(pz),
 nextParticle(NULL),daughterParticle(NULL),trackID(-1),
 mass(0.),polX(0.),polY(0.),polZ(0.)
{ PDGcode = Gcode->GetPDGEncoding(); }

G4PrimaryParticle::~G4PrimaryParticle()
{
  if(nextParticle != NULL)
  { delete nextParticle; }
  if(daughterParticle != NULL)
  { delete daughterParticle; }
}

void G4PrimaryParticle::SetPDGcode(G4int Pcode)
{
  PDGcode = Pcode;
  G4code = G4ParticleTable::GetParticleTable()->FindParticle(Pcode);
}

void G4PrimaryParticle::SetG4code(G4ParticleDefinition* Gcode)
{
  G4code = Gcode;
  PDGcode = Gcode->GetPDGEncoding();
}

const G4PrimaryParticle & 
G4PrimaryParticle::operator=(const G4PrimaryParticle &right)
{ return *this; }
int G4PrimaryParticle::operator==(const G4PrimaryParticle &right) const
{ return false; }
int G4PrimaryParticle::operator!=(const G4PrimaryParticle &right) const
{ return true; }

void G4PrimaryParticle::Print() const
{
  G4cout << "==== PDGcode " << PDGcode << "  Particle name ";
  if(G4code != NULL)
  { G4cout << G4code->GetParticleName() << endl; }
  else
  { G4cout << "is not defined in G4." << endl; }
  G4cout << "     Momentum ( " << Px << ", " << Py << ", " << Pz << " )" << endl;
  G4cout << "     Polarization ( " << polX << ", " << polY << ", "
                                 << polZ << " )" << endl;
  if(daughterParticle != NULL)
  {
    G4cout << ">>>> Daughters" << endl;
    daughterParticle->Print();
  }
  if(nextParticle != NULL)
  { nextParticle->Print(); }
  else
  { G4cout << "<<<< End of link" << endl; }
}

