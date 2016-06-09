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
#ifndef G4CASCAD_PARTICLE_HH
#define G4CASCAD_PARTICLE_HH

#include "G4InuclElementaryParticle.hh"

class G4CascadParticle {

public:

  G4CascadParticle();

  G4CascadParticle(const G4InuclElementaryParticle& particle, 
		   const std::vector<G4double>& pos,
		   G4int izone, 
		   G4double cpath,
                   G4int gen) 

    : theParticle(particle), 
    position(pos), 
    current_zone(izone), 
    current_path(cpath) {
    current_path = cpath; 
    movingIn = true;
    reflectionCounter = 0;
    generation = gen;
  };

  void updateParticleMomentum(const G4CascadeMomentum& mom) {
    theParticle.setMomentum(mom);
  };

  void updatePosition(const std::vector<G4double>& pos) {
    position = pos;
  };

  void incrementReflectionCounter() {
    reflectionCounter++; 
    reflected = true; 
  };

  void resetReflection() { 
    reflected = false; 
  };

  void incrementCurrentPath(G4double npath) { 
    current_path += npath; 
  };

  void updateZone(G4int izone) {
    current_zone = izone; 
  };

  G4bool movingInsideNuclei() const { 
    return movingIn; 
  };

  G4double getPathToTheNextZone(G4double rz_in, 
				G4double rz_out);

  const G4CascadeMomentum& getMomentum() const { 
    return theParticle.getMomentum(); 
  };

  G4InuclElementaryParticle getParticle() const { 
    return theParticle; 
  };

  const std::vector<G4double>& getPosition() const { 
    return position; 
  };

  G4int getCurrentZone() const { 
    return current_zone; 
  };

  G4int getNumberOfReflections() const { 
    return reflectionCounter; 
  };

  G4bool young(G4double young_path_cut, 
	       G4double cpath) const { 
   
    if(current_path < 1000.0) {
      return cpath < young_path_cut;
    }
    else {
      return false;
    };    
    // return current_path + cpath < young_path_cut; 
  };

  G4bool reflectedNow() const { 
    return reflected; 
  };

  void propagateAlongThePath(G4double path); 

  void print() const {
    theParticle.printParticle();
    G4cout << " zone " << current_zone << " current_path " << current_path
	   << " reflectionCounter " << reflectionCounter << G4endl
	   << " x " << position[0] << " y " << position[1]
	   << " z " << position[2] << G4endl;
  };

  G4int getGeneration() {
    return generation;
  }
   
private: 

  G4int verboseLevel;
  G4InuclElementaryParticle theParticle;
  std::vector<G4double> position;
  G4int current_zone;
  G4double current_path;
  G4bool movingIn;
  G4int reflectionCounter;   
  G4bool reflected;
  G4int generation;
};        

#endif // G4CASCAD_PARTICLE_HH
