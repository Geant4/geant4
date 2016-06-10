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
//
// $Id: G4GMocrenTouchable.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
//
// Created:  Mar. 31, 2009  Akinori Kimura  
//
// G4GMocrenTouchable class
// is used to get densities of each box cell of patient geometry
// defined by using the nest parameterisation.
//
#ifndef G4GMOCRENTOUCHABLE_HH
#define G4GMOCRENTOUCHABLE_HH

#include "G4VTouchable.hh"

class G4GMocrenTouchable : public G4VTouchable {

 public:  // with description

  G4GMocrenTouchable() {;}
  G4GMocrenTouchable(G4int & _depth0, G4int & _depth1);
  virtual ~G4GMocrenTouchable() {;}
    // Constructor and destructor.

  virtual const G4ThreeVector& GetTranslation(G4int depth=0) const;
  virtual const G4RotationMatrix* GetRotation(G4int depth=0) const;
    // Accessors for translation and rotation.

  virtual G4int GetReplicaNumber(G4int depth=0) const;
    // Methods for touchables with history.

  void SetReplicaNumber(G4int _depth0, G4int _depth1);

private:
  G4int repno[2];
};

inline
G4GMocrenTouchable::G4GMocrenTouchable(G4int & _depth0, G4int & _depth1) {
  repno[0] = _depth0;
  repno[1] = _depth1;
}

const G4ThreeVector& G4GMocrenTouchable::GetTranslation(G4int depth) const {
  // never used
  // in the purpose to avoid a warning in the compile process
  G4ThreeVector * vec = new G4ThreeVector();
  *vec *= static_cast<G4double>(depth);
  return *vec;
}
const G4RotationMatrix* G4GMocrenTouchable::GetRotation(G4int depth) const {
  // never used
  // in the puspose to avoid a warning in the compile process
  G4RotationMatrix * rot = new G4RotationMatrix();
  rot->setPhi(static_cast<G4double>(depth));
  return rot;
}
inline
G4int G4GMocrenTouchable::GetReplicaNumber(G4int depth) const {
  if(depth > 1) {
    G4Exception("G4GMocrenTouchable::GetReplicaNumber(G4int)", "gMocren0001",
		FatalException, "depth number is less than 2.");
  }
  G4int rvalue;
  if(depth < 2) rvalue = depth;
  else rvalue = 0;
  return  rvalue;
}


#endif
