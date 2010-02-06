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
#ifndef G4INTERACTION_CASE_HH
#define G4INTERACTION_CASE_HH

//#ifndef G4INUCL_PARTICLE_HH
#include "G4InuclParticle.hh"
//#endif

#include <algorithm>

class G4InteractionCase {

public:

  G4InteractionCase() { 
    bultag = std::pair<G4InuclParticle*, G4InuclParticle*>(0, 0);
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
    bultag = std::pair<G4InuclParticle*, G4InuclParticle*>(part1, part2);
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

  std::pair<G4InuclParticle*, G4InuclParticle*> bultag;

  G4int inter_case;

};

#endif // G4INTERACTION_CASE_HH 


