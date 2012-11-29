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
#ifndef G4VCASCADE_DEEXCITATION_HH
#define G4VCASCADE_DEEXCITATION_HH
// $Id$
//
// Base class to define a common interface for post-cascade processing.

#include "G4CascadeColliderBase.hh"
#include "globals.hh"
#include "G4CollisionOutput.hh"

class G4InuclParticle;
class G4Fragment;


class G4VCascadeDeexcitation : public G4CascadeColliderBase {
public:
  G4VCascadeDeexcitation(const char* name) : G4CascadeColliderBase(name) {}
  virtual ~G4VCascadeDeexcitation() {}

  // Standard Collider interface (bullet is not used in this case)
  virtual void collide(G4InuclParticle* /*bullet*/, G4InuclParticle* target,
		       G4CollisionOutput& globalOutput) = 0;

  // Interface specific to pre-compound (post-cascade) processing
  virtual void deExcite(G4Fragment* fragment,
			G4CollisionOutput& globalOutput) = 0;

protected:
  G4CollisionOutput output;	// Local buffer for de-excitation stages
};

#endif	/* G4CASCADE_DEEXCITATION_HH */
