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
#ifndef G4INTRA_NUCLEI_CASCADER_HH
#define G4INTRA_NUCLEI_CASCADER_HH

#include "G4Collider.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4CascadSpecialFunctions.hh"
#include "G4InuclElementaryParticle.hh"

using namespace G4InuclSpecialFunctions;
using namespace G4CascadSpecialFunctions;

class G4IntraNucleiCascader : public G4Collider {

public:

  G4IntraNucleiCascader();

  void setElementaryParticleCollider(G4ElementaryParticleCollider* ecollider) {
    theElementaryParticleCollider = ecollider;   
  };
  
  virtual G4CollisionOutput collide(G4InuclParticle* bullet,
				  G4InuclParticle* target);

  void setInteractionCase(G4int intcase) { 
    inter_case = intcase; 
  };

private: 
G4int verboseLevel;
  G4ElementaryParticleCollider* theElementaryParticleCollider;

  G4int inter_case;

  G4bool goodCase(G4double a, G4double z, G4double eexs, G4double ein) const; 

};        

#endif // G4INTRA_NUCLEI_CASCADER_HH 
