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
#ifndef G4V_CASCADE_COLLIDER_HH
#define G4V_CASCADE_COLLIDER_HH
// $Id: G4VCascadeCollider.hh 71719 2013-06-21 00:01:54Z mkelsey $
//
// 20100615  M. Kelsey -- Split constructor to have verbose separately
// 20100711  M. Kelsey -- Allow name to be changed after ctor, by self
// 20100714  M. Kelsey -- Move concrete functions to G4CascadeColliderBase
// 20130620  Address Coverity complaint about missing copy actions
// 20140930  Change name from "const char*" to "const G4String"

#include "globals.hh"

class G4InuclParticle;
class G4CollisionOutput;

class G4VCascadeCollider {
public:
  G4VCascadeCollider(const G4String& name, G4int verbose=0);

  virtual ~G4VCascadeCollider() {}

  virtual void collide(G4InuclParticle* bullet, G4InuclParticle* target,
		       G4CollisionOutput& output) = 0;

  virtual void setVerboseLevel(G4int verbose=0) { verboseLevel=verbose; }

protected:
  G4String theName;
  G4int verboseLevel;

  virtual void setName(const G4String& name) { theName = name; }

private:
  // Copying of modules is forbidden
  G4VCascadeCollider(const G4VCascadeCollider&);
  G4VCascadeCollider& operator=(const G4VCascadeCollider&);
};        

#endif	/* G4V_CASCADE_COLLIDER_HH */
