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
// $Id: G4PrimaryParticle.cc,v 1.7 2010-08-11 17:14:02 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4PrimaryParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ios.hh"
#include "G4VUserPrimaryParticleInformation.hh"

G4Allocator<G4PrimaryParticle> aPrimaryParticleAllocator;

G4PrimaryParticle::G4PrimaryParticle()
:PDGcode(0),G4code(0),
 direction(0.,0.,1.),kinE(0.),
 nextParticle(0),daughterParticle(0),trackID(-1),
 mass(-1.),charge(DBL_MAX),polX(0.),polY(0.),polZ(0.),
 Weight0(1.0),properTime(0.0),userInfo(0)
{;}

G4PrimaryParticle::G4PrimaryParticle(G4int Pcode)
:PDGcode(Pcode),
 direction(0.,0.,1.),kinE(0.),
 nextParticle(0),daughterParticle(0),trackID(-1),
 mass(-1.),charge(DBL_MAX),polX(0.),polY(0.),polZ(0.),
 Weight0(1.0),properTime(0.0),userInfo(0)
{ 
  G4code = G4ParticleTable::GetParticleTable()->FindParticle(Pcode); 
  if (G4code !=0) {
    mass = G4code->GetPDGMass();
  } 
}

G4PrimaryParticle::G4PrimaryParticle(G4int Pcode,
                        G4double px,G4double py,G4double pz)
:PDGcode(Pcode),
 direction(0.,0.,1.),kinE(0.),
 nextParticle(0),daughterParticle(0),trackID(-1),
 mass(-1.),charge(DBL_MAX),polX(0.),polY(0.),polZ(0.),
 Weight0(1.0),properTime(0.0),userInfo(0)
{ 
  G4code = G4ParticleTable::GetParticleTable()->FindParticle(Pcode); 
  if (G4code !=0) {
    mass = G4code->GetPDGMass();
  } 
  SetMomentum( px, py, pz);
}

G4PrimaryParticle::G4PrimaryParticle(G4int Pcode,
                        G4double px,G4double py,G4double pz,G4double E)
:PDGcode(Pcode),
 direction(0.,0.,1.),kinE(0.),
 nextParticle(0),daughterParticle(0),trackID(-1),
 charge(DBL_MAX),polX(0.),polY(0.),polZ(0.),
 Weight0(1.0),properTime(0.0),userInfo(0)
{
 G4code = G4ParticleTable::GetParticleTable()->FindParticle(Pcode); 
 Set4Momentum( px, py, pz, E); 
}

G4PrimaryParticle::G4PrimaryParticle(const G4ParticleDefinition* Gcode)
:G4code(Gcode),
 direction(0.,0.,1.),kinE(0.),
 nextParticle(0),daughterParticle(0),trackID(-1),
 mass(-1.),charge(DBL_MAX),polX(0.),polY(0.),polZ(0.),
 Weight0(1.0),properTime(0.0),userInfo(0)
{ 
  if (G4code !=0) {
    PDGcode = G4code->GetPDGEncoding(); 
    mass = G4code->GetPDGMass();
  } 
}

G4PrimaryParticle::G4PrimaryParticle(const G4ParticleDefinition* Gcode,
                        G4double px,G4double py,G4double pz)
:G4code(Gcode),
 direction(0.,0.,1.),kinE(0.),
 nextParticle(0),daughterParticle(0),trackID(-1),
 mass(-1.),charge(DBL_MAX),polX(0.),polY(0.),polZ(0.),
 Weight0(1.0),properTime(0.0),userInfo(0)
{ 
  if (G4code !=0) {
    PDGcode = G4code->GetPDGEncoding(); 
    mass = G4code->GetPDGMass();
  } 
  SetMomentum( px, py, pz);
}

G4PrimaryParticle::G4PrimaryParticle(const G4ParticleDefinition* Gcode,
                        G4double px,G4double py,G4double pz,G4double E)
:G4code(Gcode),
 direction(0.,0.,1.),kinE(0.),
 nextParticle(0),daughterParticle(0),trackID(-1),
 charge(DBL_MAX),polX(0.),polY(0.),polZ(0.),
 Weight0(1.0),properTime(0.0),userInfo(0)
{
  if (G4code !=0) {
    PDGcode = G4code->GetPDGEncoding(); 
    mass = G4code->GetPDGMass();
  } 
 Set4Momentum( px, py, pz, E); 
}

G4PrimaryParticle::~G4PrimaryParticle()
{
  if(nextParticle != 0)
  { delete nextParticle; }
  if(daughterParticle != 0)
  { delete daughterParticle; }
  if(userInfo!=0)
  { delete userInfo; }
}

void G4PrimaryParticle::SetMomentum(G4double px, G4double py, G4double pz)
{ 
  if ((mass<0.)&&(G4code!=0)){ 
    mass =  G4code->GetPDGMass(); 
  }
  G4double pmom =  sqrt(px*px+py*py+pz*pz);
  if (pmom>0.0) {
    direction.setX(px/pmom);
    direction.setY(py/pmom);
    direction.setZ(pz/pmom);
  }
  kinE = sqrt(px*px+py*py+pz*pz+mass*mass)-mass;
}

void G4PrimaryParticle::Set4Momentum(G4double px, G4double py, G4double pz, G4double E)
{ 
  G4double pmom =  sqrt(px*px+py*py+pz*pz);
  if (pmom>0.0) {
    direction.setX(px/pmom);
    direction.setY(py/pmom);
    direction.setZ(pz/pmom);
  }
  G4double mas2 = E*E - pmom*pmom;
  if(mas2>=0.){ 
    mass = std::sqrt(mas2); 
  } else { 
    if (G4code!=0){ 
      mass =  G4code->GetPDGMass(); 
    }
    E = sqrt(pmom*pmom+mass*mass);
  }
  kinE = E - mass;
}

void G4PrimaryParticle::SetPDGcode(G4int Pcode)
{
  PDGcode = Pcode;
  G4code = G4ParticleTable::GetParticleTable()->FindParticle(Pcode);
  if (G4code!=0){ 
    mass =  G4code->GetPDGMass(); 
  }
}

void G4PrimaryParticle::SetG4code(const G4ParticleDefinition* Gcode)
{
  SetParticleDefinition(Gcode);
}

void G4PrimaryParticle::SetParticleDefinition(const G4ParticleDefinition* Gcode)
{
  G4code = Gcode;
  if (G4code!=0){ 
    PDGcode = G4code->GetPDGEncoding();
    mass =  G4code->GetPDGMass(); 
  }
}

G4double G4PrimaryParticle::GetMass() const
{
  if (mass <0.0 ) {
    if (G4code != 0) {
      // return PDG mass if dynamical mass has not be specified 
      return G4code->GetPDGMass();
    }
  }
  return mass; 
}

G4double G4PrimaryParticle::GetCharge() const
{
  if ( charge <DBL_MAX ) {
    return charge;
  } else {
    if (G4code != 0) {
      // return PDG charge if dynamical mass has not be specified 
      return G4code->GetPDGCharge();
    }
  }
  return charge; 
}

const G4PrimaryParticle & 
G4PrimaryParticle::operator=(const G4PrimaryParticle &)
{ return *this; }
G4int G4PrimaryParticle::operator==(const G4PrimaryParticle &right) const
{ return (this==&right); }
G4int G4PrimaryParticle::operator!=(const G4PrimaryParticle &right) const
{ return (this!=&right); }

void G4PrimaryParticle::Print() const
{
  G4cout << "==== PDGcode " << PDGcode << "  Particle name ";
  if(G4code != 0)
  { G4cout << G4code->GetParticleName() << G4endl; }
  else
  { G4cout << " is not defined in G4." << G4endl; }
  if(charge<DBL_MAX)
  { G4cout << " Assigned charge : " << charge/eplus  << G4endl; }
  else
  { G4cout << " Charge will be taken from PDG charge." << G4endl; }
  G4cout << "     Momentum ( " 
	 << GetTotalMomentum()*direction.x()/GeV << "[GeV/c], " 
	 << GetTotalMomentum()*direction.y()/GeV << "[GeV/c], " 
	 << GetTotalMomentum()*direction.z()/GeV << "[GeV/c] )" << G4endl;
  G4cout << "     kinetic Energy : " << kinE/GeV  << " [GeV]" << G4endl;
  if(mass>=0.)
    { G4cout << "     Mass : " << mass/GeV << " [GeV]" << G4endl; }
  else
    { G4cout << "     Nominal mass will be taken from particle definition." << G4endl; }
  G4cout << "     Polarization ( " 
	 << polX << ", " 
	 << polY << ", "
	 << polZ << " )" 
	 << G4endl;
  G4cout << "     Weight : " << Weight0 << G4endl;
  if(properTime>0.0)
  { G4cout << "     PreAssigned proper decay time : " << properTime/ns << " [ns] " << G4endl; }
  if(userInfo != 0) userInfo->Print();
  if(daughterParticle != 0)
  {
    G4cout << ">>>> Daughters" << G4endl;
    daughterParticle->Print();
  }
  if(nextParticle != 0)
  { nextParticle->Print(); }
  else
  { G4cout << "<<<< End of link" << G4endl; }
}




