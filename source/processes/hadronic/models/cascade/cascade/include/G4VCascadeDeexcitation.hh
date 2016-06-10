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
// $Id: G4VCascadeDeexcitation.hh 71942 2013-06-28 19:08:11Z mkelsey $
//
// Base class to define a common interface for post-cascade processing.
//
// 20130620  Address Coverity complaint about missing copy actions
// 20130621  Implement collide() here to throw exception, use reference
//		instead of pointer for deExcite()
// 20130622  Promote to direct base of G4VCascadeCollider
// 20140930  Change name from "const char*" to "const G4String"

#include "G4VCascadeCollider.hh"
#include "globals.hh"
#include "G4CollisionOutput.hh"

class G4InuclParticle;
class G4Fragment;


class G4VCascadeDeexcitation : public G4VCascadeCollider {
public:
  G4VCascadeDeexcitation(const G4String& name) : G4VCascadeCollider(name) {}
  virtual ~G4VCascadeDeexcitation() {}

  // Standard Collider interface should not be used (will end job)
  virtual void collide(G4InuclParticle* bullet, G4InuclParticle* target,
		       G4CollisionOutput& globalOutput);

  // Interface specific to pre-compound (post-cascade) processing
  virtual void deExcite(const G4Fragment& fragment,
			G4CollisionOutput& output) = 0;

private:
  // Copying of modules is forbidden
  G4VCascadeDeexcitation(const G4VCascadeDeexcitation&);
  G4VCascadeDeexcitation& operator=(const G4VCascadeDeexcitation&);
};

#endif	/* G4CASCADE_DEEXCITATION_HH */
