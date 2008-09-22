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
#include "G4CascadParticle.hh"

G4CascadParticle::G4CascadParticle()
  : verboseLevel(2) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4CascadParticle::G4CascadParticle" << G4endl;
  }
}

G4double G4CascadParticle::getPathToTheNextZone(G4double rz_in, 
						G4double rz_out) {
  verboseLevel = 2;

  if (verboseLevel > 3) {
    G4cout << " >>> G4CascadParticle::getPathToTheNextZone" << G4endl;
  }

  G4double path = -1.0;
  G4double rp = 0.0;
  G4double rr = 0.0;
  G4double pp = 0.0;
  const G4CascadeMomentum& mom = theParticle.getMomentum();

  for (G4int i = 1; i < 4; i++) {
    rp += mom[i] * position[i - 1]; 
    rr += position[i - 1] * position[i - 1]; 
    pp += mom[i] * mom[i]; 
  };

  G4double ra = rr - rp * rp / pp;
  pp = std::sqrt(pp);
  G4double ds;
  G4double d2;
 
  if (current_zone == 0 || rp > 0.0) {
    d2 = rz_out * rz_out - ra;
    ds = 1.0;
    movingIn = false; 
 
  } else { 

    d2 = rz_in * rz_in - ra;

    if (d2 > 0.0) {
      ds = -1.0;
      movingIn = true;
  
    } else {

      d2 = rz_out * rz_out - ra;
      ds = 1.0;
      movingIn = false;
    };
  };

  path = ds * std::sqrt(d2) - rp / pp;

  return path;    
}

void G4CascadParticle::propagateAlongThePath(G4double path) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4CascadParticle::propagateAlongThePath" << G4endl;
  }

  const G4CascadeMomentum& mom = theParticle.getMomentum();
  G4double pmod = theParticle.getMomModule();

  for(G4int i = 0; i < 3; i++) position[i] += mom[i + 1] * path / pmod;

}
