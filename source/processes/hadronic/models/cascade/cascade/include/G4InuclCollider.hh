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
// $Id: G4InuclCollider.hh,v 1.22 2010-12-15 07:39:56 gunter Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100517  M. Kelsey -- Inherit from common base class, make other colliders
//		simple data members
// 20100620  M. Kelsey -- Move output buffers here to reduce memory churn
// 20100714  M. Kelsey -- Switch to new G4CascadeColliderBase class
// 20100720  M. Kelsey -- Make all the collders pointer members (to reducde
//		external compile dependences).
// 20100915  M. Kelsey -- Move de-excitation colliders to G4CascadeDeexcitation
// 20100922  M. Kelsey -- Add functions to select de-excitation method, change
//		"theDeexcitation" to be a base-class pointer for switching
// 20100926  M. Kelsey -- Use new intermediate base class for de-exciations.

#ifndef G4INUCL_COLLIDER_HH
#define G4INUCL_COLLIDER_HH

#include "G4CascadeColliderBase.hh"
#include "G4CollisionOutput.hh"

class G4InuclParticle;
class G4ElementaryParticleCollider;
class G4IntraNucleiCascader;
class G4VCascadeDeexcitation;


class G4InuclCollider : public G4CascadeColliderBase {
public:
  G4InuclCollider();
  virtual ~G4InuclCollider();

  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       G4CollisionOutput& globalOutput);

  // Select betweeen different post-cascade de-excitation models
  void useCascadeDeexcitation();
  void usePreCompoundDeexcitation();

private: 
  G4ElementaryParticleCollider* theElementaryParticleCollider;
  G4IntraNucleiCascader* theIntraNucleiCascader;

  G4VCascadeDeexcitation* theDeexcitation;	// User switchable!

  G4CollisionOutput output;		// Secondaries from main cascade
  G4CollisionOutput DEXoutput;		// Secondaries from de-excitation
};        

#endif /* G4INUCL_COLLIDER_HH */


