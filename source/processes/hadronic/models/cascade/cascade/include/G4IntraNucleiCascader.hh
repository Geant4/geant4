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
// $Id: G4IntraNucleiCascader.hh,v 1.13 2010-06-23 19:25:35 mkelsey Exp $
// GEANT4 tag: $Name: not supported by cvs2svn $
//
// 20100315  M. Kelsey -- Remove "using" directory and unnecessary #includes.
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100517  M. Kelsey -- Inherit from common base class, make other colliders
//		simple data members
// 20100617  M. Kelsey -- Make G4NucleiModel a data member, instead of
//		creating and deleting on every cycle.
// 20100623  M. Kelsey -- Undo change from 0617.  G4NucleiModel not reusable.

#ifndef G4INTRA_NUCLEI_CASCADER_HH
#define G4INTRA_NUCLEI_CASCADER_HH

#include "G4VCascadeCollider.hh"
#include "G4ElementaryParticleCollider.hh"

class G4CollisionOutput;
class G4InuclParticle;


class G4IntraNucleiCascader : public G4VCascadeCollider {
public:
  G4IntraNucleiCascader();
  virtual ~G4IntraNucleiCascader();

  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       G4CollisionOutput& output);

  // FIXME:  This should come from (or be determined by) G4InteractionCase
  void setInteractionCase(G4int intcase) { 
    inter_case = intcase; 
  };

private: 
  G4ElementaryParticleCollider theElementaryParticleCollider;

  // FIXME:  This should come from (or be determined by) G4InteractionCase
  G4int inter_case;

  G4bool goodCase(G4double a, G4double z, G4double eexs, G4double ein) const; 
};        

#endif /* G4INTRA_NUCLEI_CASCADER_HH */
