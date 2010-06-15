#ifndef G4V_CASCADE_COLLIDER_HH
#define G4V_CASCADE_COLLIDER_HH
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
// $Id: G4VCascadeCollider.hh,v 1.2 2010-06-15 22:47:25 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100615  M. Kelsey -- Split constructor to have verbose separately

#include "globals.hh"
#include "G4InteractionCase.hh"

class G4InuclNuclei;
class G4InuclParticle;
class G4CollisionOutput;

class G4VCascadeCollider {
public:
  explicit G4VCascadeCollider(const char* name);
  G4VCascadeCollider(const char* name, G4int verbose);

  virtual ~G4VCascadeCollider() {}

  virtual void collide(G4InuclParticle* bullet, G4InuclParticle* target,
		       G4CollisionOutput& output) = 0;

  virtual void setVerboseLevel(G4int verbose=0) { verboseLevel=verbose; }

protected:
  const char* theName;
  G4int verboseLevel;
  G4InteractionCase interCase;		// Determine bullet vs. target

  // Decide whether to use G4ElementaryParticleCollider or not
  virtual G4bool useEPCollider(G4InuclParticle* bullet, 
			       G4InuclParticle* target) const;

  // Decide whether to use G4BigBanger or not
  virtual G4bool explosion(G4InuclNuclei* target) const;

  // Decide whether to use G4IntraNuclearCascader or not
  virtual G4bool inelasticInteractionPossible(G4InuclParticle* bullet,
					      G4InuclParticle* target, 
					      G4double ekin) const;
};        

#endif	/* G4V_CASCADE_COLLIDER_HH */
