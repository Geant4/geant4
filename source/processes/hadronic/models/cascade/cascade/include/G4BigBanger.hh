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
#ifndef G4BIG_BANGER_HH
#define G4BIG_BANGER_HH

#include "G4Collider.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclSpecialFunctions.hh"


using namespace G4InuclSpecialFunctions;

class G4BigBanger : public G4Collider {

public:

  G4BigBanger();

  virtual G4CollisionOutput collide(G4InuclParticle* bullet,
				    G4InuclParticle* target);

private: 

G4int verboseLevel;
  G4std::vector<G4InuclElementaryParticle> generateBangInSCM(G4double etot, 
						      G4double a, 
						      G4double z, 
						      G4double mp,
						      G4double mn) const;

  G4std::vector<G4double> generateMomentumModules(G4double etot, 
					   G4double a, 
					   G4double z,
					   G4double mp, 
					   G4double mn) const; 

  G4double xProbability(G4double x, 
			G4int ia) const; 

  G4double maxProbability(G4double a) const;

  G4double generateX(G4int ia, 
		     G4double a, 
		     G4double promax) const; 

};        

#endif // G4BIG_BANGER_HH 











