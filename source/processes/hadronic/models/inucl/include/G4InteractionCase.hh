#ifndef G4INTERACTION_CASE_HH
#define G4INTERACTION_CASE_HH

#include "G4InuclParticle.hh"

#include "pair.h"

class G4InteractionCase {

public:

  G4InteractionCase() { 
    bultag = pair<G4InuclParticle*, G4InuclParticle*>(0, 0);
    inter_case = 0;
  };

  G4InteractionCase(G4InuclParticle* part1, 
		    G4InuclParticle* part2, 
		    G4int ic) {
    setBulletTarget(part1, part2);
    setInterCase(ic);
  }; 

  void setBulletTarget(G4InuclParticle* part1, 
		       G4InuclParticle* part2) {
    bultag = pair<G4InuclParticle*, G4InuclParticle*>(part1, part2);
  };

  void setInterCase(G4int ic) { 
    inter_case = ic; 
  };

  G4InuclParticle* getBullet() const { 
    return bultag.first; 
  };

  G4InuclParticle* getTarget() const { 
    return bultag.second; 
  };

  G4int getInterCase() const { 
    return inter_case; 
  };

private:

  pair<G4InuclParticle*, InuclParticle*> bultag;

  G4int inter_case;

};

#endif // G4INTERACTION_CASE_HH 


