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
#include "G4InuclElementaryParticle.hh"

class G4CascadParticle {

public:

  G4CascadParticle();

  G4CascadParticle(const G4InuclElementaryParticle& particle, 
		   const G4std::vector<G4double>& pos,
		   G4int izone, 
		   G4double cpath) 

    : theParticle(particle), 
    position(pos), 
    current_zone(izone), 
    current_path(cpath) {
    current_path = cpath; 
    movingIn = true;
    reflectionCounter = 0;   
  };

  void updateParticleMomentum(const G4std::vector<G4double>& mom) {
    theParticle.setMomentum(mom);
  };

  void updatePosition(const G4std::vector<G4double>& pos) {
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

  G4std::vector<G4double> getMomentum() const { 
    return theParticle.getMomentum(); 
  };

  G4InuclElementaryParticle getParticle() const { 
    return theParticle; 
  };

  G4std::vector<G4double> getPosition() const { 
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
   
private: 

  G4int verboseLevel;
  G4InuclElementaryParticle theParticle;
  G4std::vector<G4double> position;
  G4int current_zone;
  G4double current_path;
  G4bool movingIn;
  G4int reflectionCounter;   
  G4bool reflected;
 
};        

#endif // G4CASCAD_PARTICLE_HH
