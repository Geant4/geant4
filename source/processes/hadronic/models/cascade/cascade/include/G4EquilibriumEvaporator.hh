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
// $Id: G4EquilibriumEvaporator.hh,v 1.12 2010-06-18 02:57:44 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100517  M. Kelsey -- Inherit from common base class, make other colliders
//		simple data members.  Rename timeToBigBang() to override
//		base explosion().

#ifndef G4EQUILIBRIUM_EVAPORATOR_HH
#define G4EQUILIBRIUM_EVAPORATOR_HH

#include "G4VCascadeCollider.hh"
#include "G4Fissioner.hh"
#include "G4BigBanger.hh"

class G4CollisionOutput;
class G4InuclParticle;

class G4EquilibriumEvaporator : public G4VCascadeCollider {
public:
  G4EquilibriumEvaporator();
  virtual ~G4EquilibriumEvaporator();

  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       G4CollisionOutput& output);

private: 
  // Replace base class verision
  virtual G4bool explosion(G4InuclNuclei*) const { return false; }
  virtual G4bool explosion(G4double a, 
			   G4double z, 
			   G4double e) const;

  G4bool goodRemnant(G4double a, 
		     G4double z) const; 

  G4double getE0(G4double A) const; 

  G4double getPARLEVDEN(G4double A, 
			G4double Z) const; 

  G4double getQF(G4double x, 
		 G4double x2, 
		 G4double a, 
		 G4double z, 
		 G4double e) const;

  G4double getAF(G4double x, 
		 G4double a, 
		 G4double z, 
		 G4double e) const; 

  G4Fissioner theFissioner;
  G4BigBanger theBigBanger;
};        

#endif /* G4EQUILIBRIUM_EVAPORATOR_HH */


