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
// $Id$
//
// 20100112  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100126  M. Kelsey -- Replace vector<G4Double> position with G4ThreeVector,
//		move ::print() to .cc file, fix uninitialized data members
// 20100915  M. Kelsey -- Make getGeneration() const
// 20110729  M. Kelsey -- Add initializer for _all_ data members (path, gen),
//		re-organize declarations, with set/get pairs together
// 20110806  M. Kelsey -- Add fill() function to replicate ctor/op=() action
// 20110922  M. Kelsey -- Add stream argument to print(), add operator<<().
// 20120306  M. Kelsey -- Add access for cumulative path through nucleus.

#ifndef G4CASCAD_PARTICLE_HH
#define G4CASCAD_PARTICLE_HH

#include "G4InuclElementaryParticle.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include <iosfwd>


class G4CascadParticle {

public:
  // NOTE:  Default constructor does not make a functional object!
  G4CascadParticle();

  G4CascadParticle(const G4InuclElementaryParticle& particle, 
		   const G4ThreeVector& pos, G4int izone, G4double cpath,
                   G4int gen) 
    : verboseLevel(0), theParticle(particle), position(pos), 
      current_zone(izone), current_path(cpath), movingIn(true),
      reflectionCounter(0), reflected(false), generation(gen) {}

  // Analogue to operator=() to support filling vectors w/o temporaries
  void fill(const G4InuclElementaryParticle& particle, 
	    const G4ThreeVector& pos, G4int izone, G4double cpath,
	    G4int gen);

  const G4InuclElementaryParticle& getParticle() const { return theParticle; }
  G4InuclElementaryParticle& getParticle() { return theParticle; }

  G4int getGeneration() const { return generation; }
  void setGeneration(G4int gen) { generation = gen; }

  G4LorentzVector getMomentum() const {		// Can't return ref; temporary
    return theParticle.getMomentum(); 
  }

  void updateParticleMomentum(const G4LorentzVector& mom) {
    theParticle.setMomentum(mom);
  }

  const G4ThreeVector& getPosition() const { return position; }
  void updatePosition(const G4ThreeVector& pos) { position = pos; }

  void incrementReflectionCounter() {
    reflectionCounter++; 
    reflected = true; 
  }
  G4int getNumberOfReflections() const { return reflectionCounter; }

  void resetReflection() { reflected = false; }
  G4bool reflectedNow() const { return reflected; }

  void initializePath(G4double npath) { current_path = npath; }
  void incrementCurrentPath(G4double npath) { current_path += npath; }
  G4double getCurrentPath() const { return current_path; }

  void updateZone(G4int izone) { current_zone = izone; }
  G4int getCurrentZone() const { return current_zone; }

  void setMovingInsideNuclei(G4bool isMovingIn=true) { movingIn = isMovingIn; }
  G4bool movingInsideNuclei() const { return movingIn; }

  G4double getPathToTheNextZone(G4double rz_in, G4double rz_out);
  void propagateAlongThePath(G4double path);

  G4bool young(G4double young_path_cut, G4double cpath) const { 
    return ((current_path < 1000.) && (cpath < young_path_cut));
  }

  void print(std::ostream& os) const;

private: 
  G4int verboseLevel;
  G4InuclElementaryParticle theParticle;
  G4ThreeVector position;
  G4int current_zone;
  G4double current_path;
  G4bool movingIn;
  G4int reflectionCounter;   
  G4bool reflected;
  G4int generation;
};        

// Proper stream output (just calls print())

std::ostream& operator<<(std::ostream& os, const G4CascadParticle& part);

#endif // G4CASCAD_PARTICLE_HH
