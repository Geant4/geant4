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
#ifndef G4FISSIONER_HH
#define G4FISSIONER_HH

#include "G4Collider.hh"
#include "G4InuclSpecialFunctions.hh"

using namespace G4InuclSpecialFunctions;

class G4Fissioner : public G4Collider {

public:

  G4Fissioner();

  virtual G4CollisionOutput collide(G4InuclParticle* bullet,
				    G4InuclParticle* target);

private: 
G4int verboseLevel;
  G4double getC2(G4double A1, 
		 G4double A2, 
		 G4double X3, 
		 G4double X4, 
		 G4double R12) const; 

  G4double getZopt(G4double A1, 
		   G4double A2, 
		   G4double ZT, 
                   G4double X3, 
		   G4double X4, 
		   G4double R12) const;
		    
  void potentialMinimization(G4double& VP, 
			     G4std::vector<G4double>& ED, 
			     G4double& VC,
			     G4double AF, 
			     G4double AS, 
			     G4double ZF, 
			     G4double ZS,
			     G4std::vector<G4double>& AL1, 
			     G4std::vector<G4double>& BET1, 
			     G4double& R12) const; 

};        

#endif // G4FISSIONER_HH
