// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PrimaryParticle.cc,v 1.4 2001-02-07 08:20:43 asaim Exp $
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
 mass(0.),polX(0.),polY(0.),polZ(0.),Weight0(1.0),properTime(0.0)
{;}

G4PrimaryParticle::G4PrimaryParticle(G4int Pcode)
:PDGcode(Pcode),Px(0.),Py(0.),Pz(0.),
 nextParticle(NULL),daughterParticle(NULL),trackID(-1),
 mass(0.),polX(0.),polY(0.),polZ(0.),Weight0(1.0),properTime(0.0)
{ G4code = G4ParticleTable::GetParticleTable()->FindParticle(Pcode); }

G4PrimaryParticle::G4PrimaryParticle(G4int Pcode,
                        G4double px,G4double py,G4double pz)
:PDGcode(Pcode),Px(px),Py(py),Pz(pz),
 nextParticle(NULL),daughterParticle(NULL),trackID(-1),
 mass(0.),polX(0.),polY(0.),polZ(0.),Weight0(1.0),properTime(0.0)
{ G4code = G4ParticleTable::GetParticleTable()->FindParticle(Pcode); }

G4PrimaryParticle::G4PrimaryParticle(G4ParticleDefinition* Gcode)
:G4code(Gcode),Px(0.),Py(0.),Pz(0.),
 nextParticle(NULL),daughterParticle(NULL),trackID(-1),
 mass(0.),polX(0.),polY(0.),polZ(0.),Weight0(1.0),properTime(0.0)
{ PDGcode = Gcode->GetPDGEncoding(); }

G4PrimaryParticle::G4PrimaryParticle(G4ParticleDefinition* Gcode,
                        G4double px,G4double py,G4double pz)
:G4code(Gcode),Px(px),Py(py),Pz(pz),
 nextParticle(NULL),daughterParticle(NULL),trackID(-1),
 mass(0.),polX(0.),polY(0.),polZ(0.),Weight0(1.0),properTime(0.0)
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
  { G4cout << G4code->GetParticleName() << G4endl; }
  else
  { G4cout << "is not defined in G4." << G4endl; }
  G4cout << "     Momentum ( " << Px << ", " << Py << ", " << Pz << " )" << G4endl;
  G4cout << "     Polarization ( " << polX << ", " << polY << ", "
                                 << polZ << " )" << G4endl;
  G4cout << "     Weight : " << Weight0 << G4endl;
  if(properTime>0.0)
  { G4cout << "     PreAssigned proper decay time : " << properTime/ns << " (nsec)" << G4endl; }
  if(daughterParticle != NULL)
  {
    G4cout << ">>>> Daughters" << G4endl;
    daughterParticle->Print();
  }
  if(nextParticle != NULL)
  { nextParticle->Print(); }
  else
  { G4cout << "<<<< End of link" << G4endl; }
}

