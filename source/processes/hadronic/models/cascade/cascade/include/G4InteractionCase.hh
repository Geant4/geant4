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
#ifndef G4INTERACTION_CASE_HH
#define G4INTERACTION_CASE_HH

//#ifndef G4INUCL_PARTICLE_HH
#include "G4InuclParticle.hh"
//#endif

#include "g4std/algorithm"

class G4InteractionCase {

public:

  G4InteractionCase() { 
    bultag = G4std::pair<G4InuclParticle*, G4InuclParticle*>(0, 0);
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
    bultag = G4std::pair<G4InuclParticle*, G4InuclParticle*>(part1, part2);
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

  G4std::pair<G4InuclParticle*, G4InuclParticle*> bultag;

  G4int inter_case;

};

#endif // G4INTERACTION_CASE_HH 


