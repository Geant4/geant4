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
// INCL++ intra-nuclear cascade model
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/*
 * Particle.cc
 *
 *  \date Jun 5, 2009
 * \author Pekka Kaitaniemi
 */

#include "G4INCLParticle.hh"
#include "G4INCLParticleTable.hh"

namespace G4INCL {

#ifdef INCLXX_IN_GEANT4_MODE
    std::vector<G4double> Particle::INCLBiasVector;
#else
    G4ThreadLocal std::vector<G4double> Particle::INCLBiasVector;
  //G4VectorCache<G4double> Particle::INCLBiasVector;
#endif
    G4ThreadLocal long Particle::nextID = 1;
  G4ThreadLocal G4int Particle::nextBiasedCollisionID = 0;

  Particle::Particle()
    : theZ(0), theA(0), theS(0),
    theParticipantType(TargetSpectator),
    theType(UnknownParticle),
    theEnergy(0.0),
    thePropagationEnergy(&theEnergy),
    theFrozenEnergy(theEnergy),
    theMomentum(ThreeVector(0.,0.,0.)),
    thePropagationMomentum(&theMomentum),
    theFrozenMomentum(theMomentum),
    thePosition(ThreeVector(0.,0.,0.)),
    nCollisions(0),
    nDecays(0),
    thePotentialEnergy(0.0),
    rpCorrelated(false),
    uncorrelatedMomentum(0.),
    theParticleBias(1.),
    theNKaon(0),
    theParentResonancePDGCode(0),
    theParentResonanceID(0),
    theHelicity(0.0),
    emissionTime(0.0),
    outOfWell(false),
    theMass(0.)
  {
    ID = nextID;
    nextID++;
  }

  Particle::Particle(ParticleType t, G4double energy,
      ThreeVector const &momentum, ThreeVector const &position)
    : theEnergy(energy),
    thePropagationEnergy(&theEnergy),
    theFrozenEnergy(theEnergy),
    theMomentum(momentum),
    thePropagationMomentum(&theMomentum),
    theFrozenMomentum(theMomentum),
    thePosition(position),
    nCollisions(0), nDecays(0),
    thePotentialEnergy(0.),
    rpCorrelated(false),
    uncorrelatedMomentum(theMomentum.mag()),
    theParticleBias(1.),
    theNKaon(0),
    theParentResonancePDGCode(0),
    theParentResonanceID(0),
    theHelicity(0.0),
    emissionTime(0.0), outOfWell(false)
  {
    theParticipantType = TargetSpectator;
    ID = nextID;
    nextID++;
    if(theEnergy <= 0.0) {
      INCL_WARN("Particle with energy " << theEnergy << " created." << '\n');
    }
    setType(t);
    setMass(getInvariantMass());
  }

  Particle::Particle(ParticleType t,
      ThreeVector const &momentum, ThreeVector const &position)
    : thePropagationEnergy(&theEnergy),
    theMomentum(momentum),
    thePropagationMomentum(&theMomentum),
    theFrozenMomentum(theMomentum),
    thePosition(position),
    nCollisions(0), nDecays(0),
    thePotentialEnergy(0.),
    rpCorrelated(false),
    uncorrelatedMomentum(theMomentum.mag()),
    theParticleBias(1.),
    theNKaon(0),
    theParentResonancePDGCode(0),
    theParentResonanceID(0),
    theHelicity(0.0),
    emissionTime(0.0), outOfWell(false)
  {
    theParticipantType = TargetSpectator;
    ID = nextID;
    nextID++;
    setType(t);
    if( isResonance() ) {
      INCL_ERROR("Cannot create resonance without specifying its momentum four-vector." << '\n');
    }
    G4double energy = std::sqrt(theMomentum.mag2() + theMass*theMass);
    theEnergy = energy;
    theFrozenEnergy = theEnergy;
  }

  const ThreeVector &Particle::adjustMomentumFromEnergy() {
    const G4double p2 = theMomentum.mag2();
    G4double newp2 = theEnergy*theEnergy - theMass*theMass;
    if( newp2<0.0 ) {
      INCL_ERROR("Particle has E^2 < m^2." << '\n' << print());
      newp2 = 0.0;
      theEnergy = theMass;
    }

    theMomentum *= std::sqrt(newp2/p2);
    return theMomentum;
  }

  G4double Particle::adjustEnergyFromMomentum() {
    theEnergy = std::sqrt(theMomentum.mag2() + theMass*theMass);
    return theEnergy;
  }

  void ParticleList::rotatePositionAndMomentum(const G4double angle, const ThreeVector &axis) const {
    for(const_iterator i=begin(), e=end(); i!=e; ++i) {
      (*i)->rotatePositionAndMomentum(angle, axis);
    }
  }

  void ParticleList::rotatePosition(const G4double angle, const ThreeVector &axis) const {
    for(const_iterator i=begin(), e=end(); i!=e; ++i) {
      (*i)->rotatePosition(angle, axis);
    }
  }

  void ParticleList::rotateMomentum(const G4double angle, const ThreeVector &axis) const {
    for(const_iterator i=begin(), e=end(); i!=e; ++i) {
      (*i)->rotateMomentum(angle, axis);
    }
  }

  void ParticleList::boost(const ThreeVector &b) const {
    for(const_iterator i=begin(), e=end(); i!=e; ++i) {
      (*i)->boost(b);
    }
  }

  G4double ParticleList::getParticleListBias() const {
    if(G4int((*this).size())==0) return 1.;
    std::vector<G4int> MergedVector;
    for(ParticleIter i = (*this).begin(), e = (*this).end(); i!=e; ++i){
        MergedVector = Particle::MergeVectorBias(MergedVector,(*i));
    }
    return Particle::getBiasFromVector(MergedVector);
  }

  std::vector<G4int> ParticleList::getParticleListBiasVector() const {
    std::vector<G4int> MergedVector;
    if(G4int((*this).size())==0) return MergedVector;
    for(ParticleIter i = (*this).begin(), e = (*this).end(); i!=e; ++i){
        MergedVector = Particle::MergeVectorBias(MergedVector,(*i));
    }
    return MergedVector;
  }

  void Particle::FillINCLBiasVector(G4double newBias){
// assert(G4int(Particle::INCLBiasVector.size())==nextBiasedCollisionID);
    //assert(G4int(Particle::INCLBiasVector.Size())==nextBiasedCollisionID);
// assert(std::fabs(newBias - 1.) > 1E-6);
	Particle::INCLBiasVector.push_back(newBias);
	//Particle::INCLBiasVector.Push_back(newBias);
    Particle::nextBiasedCollisionID++;
  }

  G4double Particle::getBiasFromVector(std::vector<G4int> VectorBias) {
    if(VectorBias.empty()) return 1.;
    
    G4double ParticleBias = 1.;
    
    for(G4int i=0; i<G4int(VectorBias.size()); i++){
        ParticleBias *= Particle::INCLBiasVector[G4int(VectorBias[i])];
    }
    
    return ParticleBias;
  }
  
  std::vector<G4int> Particle::MergeVectorBias(Particle const * const p1, Particle const * const p2){
    std::vector<G4int> MergedVectorBias;
    std::vector<G4int> VectorBias1 = p1->getBiasCollisionVector();
    std::vector<G4int> VectorBias2 = p2->getBiasCollisionVector();
    G4int i = 0;
    G4int j = 0;
    if(VectorBias1.size()==0 && VectorBias2.size()==0) return MergedVectorBias;
    else if(VectorBias1.size()==0) return VectorBias2;
    else if(VectorBias2.size()==0) return VectorBias1;
        
    while(i < G4int(VectorBias1.size()) || j < G4int(VectorBias2.size())){
        if(VectorBias1[i]==VectorBias2[j]){
            MergedVectorBias.push_back(VectorBias1[i]);
            i++;
            j++;
            if(i == G4int(VectorBias1.size())){
                for(;j<G4int(VectorBias2.size());j++) MergedVectorBias.push_back(VectorBias2[j]);
            }
            else if(j == G4int(VectorBias2.size())){
                for(;i<G4int(VectorBias1.size());i++) MergedVectorBias.push_back(VectorBias1[i]);
            }
        } else if(VectorBias1[i]<VectorBias2[j]){
            MergedVectorBias.push_back(VectorBias1[i]);
            i++;
            if(i == G4int(VectorBias1.size())){
                for(;j<G4int(VectorBias2.size());j++) MergedVectorBias.push_back(VectorBias2[j]);
            }
        }
        else {
            MergedVectorBias.push_back(VectorBias2[j]);
            j++;
            if(j == G4int(VectorBias2.size())){
                for(;i<G4int(VectorBias1.size());i++) MergedVectorBias.push_back(VectorBias1[i]);
            }
        }
    }
    return MergedVectorBias;
  }
  
  std::vector<G4int> Particle::MergeVectorBias(std::vector<G4int> p1, Particle const * const p2){
    std::vector<G4int> MergedVectorBias;
    std::vector<G4int> VectorBias = p2->getBiasCollisionVector();
    G4int i = 0;
    G4int j = 0;
    if(p1.size()==0 && VectorBias.size()==0) return MergedVectorBias;
    else if(p1.size()==0) return VectorBias;
    else if(VectorBias.size()==0) return p1;
        
    while(i < G4int(p1.size()) || j < G4int(VectorBias.size())){
        if(p1[i]==VectorBias[j]){
            MergedVectorBias.push_back(p1[i]);
            i++;
            j++;
            if(i == G4int(p1.size())){
                for(;j<G4int(VectorBias.size());j++) MergedVectorBias.push_back(VectorBias[j]);
            }
            else if(j == G4int(VectorBias.size())){
                for(;i<G4int(p1.size());i++) MergedVectorBias.push_back(p1[i]);
            }
        } else if(p1[i]<VectorBias[j]){
            MergedVectorBias.push_back(p1[i]);
            i++;
            if(i == G4int(p1.size())){
                for(;j<G4int(VectorBias.size());j++) MergedVectorBias.push_back(VectorBias[j]);
            }
        }
        else {
            MergedVectorBias.push_back(VectorBias[j]);
            j++;
            if(j == G4int(VectorBias.size())){
                for(;i<G4int(p1.size());i++) MergedVectorBias.push_back(p1[i]);
            }
        }
    }
    return MergedVectorBias;
  }

  G4double Particle::getTotalBias() {
      G4double TotalBias = 1.;
      for(G4int i=0; i<G4int(INCLBiasVector.size());i++) TotalBias *= Particle::INCLBiasVector[i];
      return TotalBias;
  }

  void Particle::setINCLBiasVector(std::vector<G4double> NewVector) {
      Particle::INCLBiasVector = NewVector;
  }
}
