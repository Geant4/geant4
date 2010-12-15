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
// $Id: G4EquilibriumEvaporator.hh,v 1.16 2010-12-15 07:39:48 gunter Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100517  M. Kelsey -- Inherit from common base class, make other colliders
//		simple data members.  Rename timeToBigBang() to override
//		base explosion().
// 20100714  M. Kelsey -- Switch to new G4CascadeColliderBase class
// 20100923  M. Kelsey -- Migrate to integer A and Z
// 20100925  M. Kelsey -- Remove no longer necessary explosion() interface

#ifndef G4EQUILIBRIUM_EVAPORATOR_HH
#define G4EQUILIBRIUM_EVAPORATOR_HH

#include "G4CascadeColliderBase.hh"
#include "G4Fissioner.hh"
#include "G4BigBanger.hh"

class G4CollisionOutput;
class G4InuclParticle;

class G4EquilibriumEvaporator : public G4CascadeColliderBase {
public:
  G4EquilibriumEvaporator();
  virtual ~G4EquilibriumEvaporator();

  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       G4CollisionOutput& output);

private: 
  // Replace base class verision
  virtual G4bool explosion(G4int a, G4int z, G4double e) const;

  // FIXME:  Need to redeclare and call through base-class polymorphisms
  virtual G4bool explosion(G4InuclNuclei* target) const {
    return G4CascadeColliderBase::explosion(target);
  }

  virtual G4bool explosion(G4Fragment* target) const {
    return G4CascadeColliderBase::explosion(target);
  }

  G4bool goodRemnant(G4int a, G4int z) const; 
  G4double getE0(G4int A) const; 
  G4double getPARLEVDEN(G4int A, G4int Z) const; 
  G4double getQF(G4double x, G4double x2, G4int a, G4int z, G4double e) const;
  G4double getAF(G4double x, G4int a, G4int z, G4double e) const; 

  G4Fissioner theFissioner;
  G4BigBanger theBigBanger;
};        

#endif /* G4EQUILIBRIUM_EVAPORATOR_HH */


