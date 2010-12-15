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
// $Id: G4IntraNucleiCascader.hh,v 1.22 2010-12-15 07:39:54 gunter Exp $
// GEANT4 tag: $Name: not supported by cvs2svn $
//
// 20100315  M. Kelsey -- Remove "using" directory and unnecessary #includes.
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100517  M. Kelsey -- Inherit from common base class, make other colliders
//		simple data members
// 20100617  M. Kelsey -- Make G4NucleiModel a data member, instead of
//		creating and deleting on every cycle.
// 20100623  M. Kelsey -- Undo change from 0617.  G4NucleiModel not reusable.
// 20100714  M. Kelsey -- Switch to new G4CascadeColliderBase class
// 20100716  M. Kelsey -- Eliminate inter_case; use base-class functionality,
//		add function to compute recoil nuclear mass on the fly
// 20100720  M. Kelsey -- Make EPCollider pointer member
// 20100722  M. Kelsey -- Move cascade output buffers to .hh file
// 20100728  M. Kelsey -- Move G4NucleiModel here, as pointer member
// 20100907  M. Kelsey -- Add new "makeResidualFragment" to create
//		G4InuclNuclei at current stage of cascade
// 20100909  M. Kelsey -- Drop makeResidualFragment(), getResidualMass() and
//		local G4InuclNuclei object, replace with new RecoilMaker.
//		Move goodCase() to RecoilMaker.
// 20100916  M. Kelsey -- Add functions to handle trapped particles, and to
//		decay hyperons.

#ifndef G4INTRA_NUCLEI_CASCADER_HH
#define G4INTRA_NUCLEI_CASCADER_HH

#include "G4CascadeColliderBase.hh"
#include "G4CollisionOutput.hh"
#include <vector>

class G4CascadParticle;
class G4CascadeRecoilMaker;
class G4CollisionOutput;
class G4ElementaryParticleCollider;
class G4InuclElementaryParticle;
class G4InuclParticle;
class G4NucleiModel;


class G4IntraNucleiCascader : public G4CascadeColliderBase {
public:
  G4IntraNucleiCascader();
  virtual ~G4IntraNucleiCascader();

  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       G4CollisionOutput& output);

protected:
  void processTrappedParticle(const G4CascadParticle& trapped);
  void decayTrappedParticle(const G4CascadParticle& trapped);

private: 
  G4NucleiModel* model;
  G4ElementaryParticleCollider* theElementaryParticleCollider;
  G4CascadeRecoilMaker* theRecoilMaker;

  // Buffers for collecting result of cascade (reset on each iteration)
  G4CollisionOutput output;
  std::vector<G4CascadParticle> cascad_particles;
  std::vector<G4CascadParticle> new_cascad_particles;
  std::vector<G4InuclElementaryParticle> output_particles;
  G4ExitonConfiguration theExitonConfiguration;
};        

#endif /* G4INTRA_NUCLEI_CASCADER_HH */
