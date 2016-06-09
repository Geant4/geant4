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
// $Id$
//
// 20100315  M. Kelsey -- Remove "using" directive and unnecessary #includes.
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100517  M. Kelsey -- Inherit from common base class
// 20100714  M. Kelsey -- Switch to new G4CascadeColliderBase class
// 20100914  M. Kelsey -- Migrate to integer A and Z

#ifndef G4NON_EQUILIBRIUM_EVAPORATOR_HH
#define G4NON_EQUILIBRIUM_EVAPORATOR_HH

#include "G4CascadeColliderBase.hh"

class G4CollisionOutput;
class G4InuclParticle;

class G4NonEquilibriumEvaporator : public G4CascadeColliderBase {
public:
  G4NonEquilibriumEvaporator();
  virtual ~G4NonEquilibriumEvaporator() {}

  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       G4CollisionOutput& output);

private: 
  G4double getMatrixElement(G4int A) const;
  G4double getE0(G4int A) const; 
  G4double getParLev(G4int A, G4int Z) const;
};

#endif /* G4NON_EQUILIBRIUM_EVAPORATOR_HH */
