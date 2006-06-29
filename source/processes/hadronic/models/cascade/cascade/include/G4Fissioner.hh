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
			     std::vector<G4double>& ED, 
			     G4double& VC,
			     G4double AF, 
			     G4double AS, 
			     G4double ZF, 
			     G4double ZS,
			     std::vector<G4double>& AL1, 
			     std::vector<G4double>& BET1, 
			     G4double& R12) const; 

};        

#endif // G4FISSIONER_HH
