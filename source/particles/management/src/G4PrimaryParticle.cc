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
// $Id: G4PrimaryParticle.cc 81373 2014-05-27 13:06:32Z gcosmo $
//

#include "G4PrimaryParticle.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ios.hh"
#include "G4VUserPrimaryParticleInformation.hh"

G4ThreadLocal G4Allocator<G4PrimaryParticle> *aPrimaryParticleAllocator = 0;

G4PrimaryParticle::G4PrimaryParticle()
:PDGcode(0),G4code(0),
 direction(0.,0.,1.),kinE(0.),
 nextParticle(0),daughterParticle(0),trackID(-1),
 mass(-1.),charge(0.),polX(0.),polY(0.),polZ(0.),
 Weight0(1.0),properTime(0.0),userInfo(0)
{;}

G4PrimaryParticle::G4PrimaryParticle(G4int Pcode)
:PDGcode(Pcode),
 direction(0.,0.,1.),kinE(0.),
 nextParticle(0),daughterParticle(0),trackID(-1),
 mass(-1.),charge(0.),polX(0.),polY(0.),polZ(0.),
 Weight0(1.0),properTime(0.0),userInfo(0)
{ 
  G4code = G4ParticleTable::GetParticleTable()->FindParticle(Pcode); 
  if (G4code !=0) {
    mass   = G4code->GetPDGMass();
    charge = G4code->GetPDGCharge();
  } 
}

G4PrimaryParticle::G4PrimaryParticle(G4int Pcode,
                        G4double px,G4double py,G4double pz)
:PDGcode(Pcode),
 direction(0.,0.,1.),kinE(0.),
 nextParticle(0),daughterParticle(0),trackID(-1),
 mass(-1.),charge(0.),polX(0.),polY(0.),polZ(0.),
 Weight0(1.0),properTime(0.0),userInfo(0)
{ 
  G4code = G4ParticleTable::GetParticleTable()->FindParticle(Pcode); 
  if (G4code !=0) {
    mass   = G4code->GetPDGMass();
    charge = G4code->GetPDGCharge();
  } 
  SetMomentum( px, py, pz);
}

G4PrimaryParticle::G4PrimaryParticle(G4int Pcode,
                        G4double px,G4double py,G4double pz,G4double E)
:PDGcode(Pcode),
 direction(0.,0.,1.),kinE(0.),
 nextParticle(0),daughterParticle(0),trackID(-1),
 mass(-1.),charge(0.),polX(0.),polY(0.),polZ(0.),
 Weight0(1.0),properTime(0.0),userInfo(0)
{
 G4code = G4ParticleTable::GetParticleTable()->FindParticle(Pcode); 
 if (G4code !=0) {
    mass = G4code->GetPDGMass();
    charge = G4code->GetPDGCharge();
  } 
 Set4Momentum( px, py, pz, E); 
}

G4PrimaryParticle::G4PrimaryParticle(const G4ParticleDefinition* Gcode)
:PDGcode(0),G4code(Gcode),
 direction(0.,0.,1.),kinE(0.),
 nextParticle(0),daughterParticle(0),trackID(-1),
 mass(-1.),charge(0.),polX(0.),polY(0.),polZ(0.),
 Weight0(1.0),properTime(0.0),userInfo(0)
{ 
  if (G4code !=0) {
    PDGcode = Gcode->GetPDGEncoding(); 
    mass = G4code->GetPDGMass();
    charge = G4code->GetPDGCharge();
  } 
}

G4PrimaryParticle::G4PrimaryParticle(const G4ParticleDefinition* Gcode,
                        G4double px,G4double py,G4double pz)
:PDGcode(0),G4code(Gcode),
 direction(0.,0.,1.),kinE(0.),
 nextParticle(0),daughterParticle(0),trackID(-1),
 mass(-1.),charge(0.),polX(0.),polY(0.),polZ(0.),
 Weight0(1.0),properTime(0.0),userInfo(0)
{ 
  if (G4code !=0) {
    PDGcode = Gcode->GetPDGEncoding(); 
    mass = G4code->GetPDGMass();
    charge = G4code->GetPDGCharge();
  } 
  SetMomentum( px, py, pz);
}

G4PrimaryParticle::G4PrimaryParticle(const G4ParticleDefinition* Gcode,
                        G4double px,G4double py,G4double pz,G4double E)
:PDGcode(0),G4code(Gcode),
 direction(0.,0.,1.),kinE(0.),
 nextParticle(0),daughterParticle(0),trackID(-1),
 mass(-1.),charge(0.),polX(0.),polY(0.),polZ(0.),
 Weight0(1.0),properTime(0.0),userInfo(0)
{
  if (G4code !=0) {
    PDGcode = Gcode->GetPDGEncoding(); 
    mass = G4code->GetPDGMass();
    charge = G4code->GetPDGCharge();
  } 
  Set4Momentum( px, py, pz, E); 
}

G4PrimaryParticle::G4PrimaryParticle(const G4PrimaryParticle& right)
:PDGcode(0),G4code(0),
 direction(0.,0.,1.),kinE(0.),
 nextParticle(0),daughterParticle(0),trackID(-1),
 mass(-1.),charge(0.),polX(0.),polY(0.),polZ(0.),
 Weight0(1.0),properTime(0.0),userInfo(0)
{
  *this = right;
}

G4PrimaryParticle & G4PrimaryParticle::operator=(const G4PrimaryParticle & right)
{ 
  if (this != &right) {
    PDGcode      = right.PDGcode;
    G4code       = right.G4code;
    direction    = right.direction;
    kinE         = right.kinE;
    if (nextParticle !=0) delete nextParticle;
    if ( right.nextParticle ==0 ){
      nextParticle = 0;
    } else {
      nextParticle = new G4PrimaryParticle(*right.nextParticle);
    }
    if (daughterParticle !=0) delete daughterParticle;
    if ( right.daughterParticle ==0 ){
      daughterParticle = 0;
    } else {
      daughterParticle = new G4PrimaryParticle(*right.daughterParticle);
    }
    trackID      = right.trackID;
    mass         = right.mass;
    charge       = right.charge;
    polX         = right.polX;
    polY         = right.polY;
    polZ         = right.polZ;
    Weight0      = right.Weight0;
    properTime   = right.properTime;

    // userInfo can not be copied
    userInfo = 0;
  }
  
  return *this; 
}

G4int G4PrimaryParticle::operator==(const G4PrimaryParticle &right) const
{ return (this==&right); }

G4int G4PrimaryParticle::operator!=(const G4PrimaryParticle &right) const
{ return (this!=&right); }

G4PrimaryParticle::~G4PrimaryParticle()
{
  if(nextParticle != 0){ 
    delete nextParticle;
    nextParticle = 0;
  }
  if(daughterParticle != 0){
    delete daughterParticle; 
    daughterParticle =0;
  }
  if(userInfo!=0) {
    delete userInfo; 
    userInfo=0;
  }
}

void G4PrimaryParticle::SetMomentum(G4double px, G4double py, G4double pz)
{ 
  if ((mass<0.)&&(G4code!=0)){ 
    mass =  G4code->GetPDGMass(); 
  }
  G4double pmom =  std::sqrt(px*px+py*py+pz*pz);
  if (pmom>0.0) {
    direction.setX(px/pmom);
    direction.setY(py/pmom);
    direction.setZ(pz/pmom);
  }
  kinE = std::sqrt(px*px+py*py+pz*pz+mass*mass)-mass;
}

void G4PrimaryParticle::Set4Momentum(G4double px, G4double py, G4double pz, G4double E)
{ 
  G4double pmom =  std::sqrt(px*px+py*py+pz*pz);
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
    E = std::sqrt(pmom*pmom+mass*mass);
  }
  kinE = E - mass;
}

void G4PrimaryParticle::SetPDGcode(G4int Pcode)
{
  PDGcode = Pcode;
  G4code = G4ParticleTable::GetParticleTable()->FindParticle(Pcode);
  if (G4code!=0){ 
    mass =  G4code->GetPDGMass(); 
    charge = G4code->GetPDGCharge();
  }
}

void G4PrimaryParticle::SetParticleDefinition(const G4ParticleDefinition* Gcode)
{
  G4code = Gcode;
  if (G4code!=0){ 
    PDGcode = Gcode->GetPDGEncoding();
    mass =  G4code->GetPDGMass(); 
    charge = G4code->GetPDGCharge();
  }
}

void G4PrimaryParticle::Print() const
{
  G4cout << "==== PDGcode " << PDGcode << "  Particle name ";
  if(G4code != 0)
  { G4cout << G4code->GetParticleName() << G4endl; }
  else
  { G4cout << " is not defined in G4." << G4endl; }
  G4cout << " Assigned charge : " << charge/eplus  << G4endl; 
  G4cout << "     Momentum ( " 
	 << GetTotalMomentum()*direction.x()/GeV << "[GeV/c], " 
	 << GetTotalMomentum()*direction.y()/GeV << "[GeV/c], " 
	 << GetTotalMomentum()*direction.z()/GeV << "[GeV/c] )" << G4endl;
  G4cout << "     kinetic Energy : " << kinE/GeV  << " [GeV]" << G4endl;
  if(mass>=0.){ 
    G4cout << "     Mass : " << mass/GeV << " [GeV]" << G4endl; 
  } else { 
    G4cout << "     Mass is not assigned " << G4endl; 
  }
  G4cout << "     Polarization ( " 
	 << polX << ", " 
	 << polY << ", "
	 << polZ << " )" 
	 << G4endl;
  G4cout << "     Weight : " << Weight0 << G4endl;
  if(properTime>0.0) { 
    G4cout << "     PreAssigned proper decay time : " << properTime/ns << " [ns] " << G4endl; 
  }
  if(userInfo != 0) { userInfo->Print(); }
  if(daughterParticle != 0) {
    G4cout << ">>>> Daughters" << G4endl;
    daughterParticle->Print();
  }
  if(nextParticle != 0) { 
    nextParticle->Print(); 
  } else { 
    G4cout << "<<<< End of link" << G4endl; 
  }
}
