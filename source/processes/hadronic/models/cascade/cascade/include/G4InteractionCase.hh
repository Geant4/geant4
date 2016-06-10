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
// $Id: G4InteractionCase.hh 66241 2012-12-13 18:34:42Z gunter $
//
// 20100518  M. Kelsey -- Why use std::pair<> at all?  Never exported; just
//		store pointers.  Add clear() function.  Move code from
//		Colliders' "bulletTargetSetter()" to set().

#ifndef G4INTERACTION_CASE_HH
#define G4INTERACTION_CASE_HH

#include "globals.hh"

class G4InuclParticle;


class G4InteractionCase {
public:
  G4InteractionCase() : bullet(0), target(0), inter_case(0) {}

  G4InteractionCase(G4InuclParticle* part1, G4InuclParticle* part2) {
    set(part1, part2);
  }

  void set(G4InuclParticle* part1, G4InuclParticle* part2);

  void clear() {
    bullet = target = 0;
    inter_case = 0;
  }

  G4InuclParticle* getBullet() const { return bullet; }
  G4InuclParticle* getTarget() const { return target; }

  G4bool valid() const      { return inter_case != 0; }

  G4bool twoNuclei() const  { return inter_case == -2; }
  G4bool hadNucleus() const { return inter_case == -1; }
  G4int  hadrons() const    { return inter_case; }	// "rtype" or "is" code

  // For compatibility with G4IntraNucleiCascader code
  G4int  code() const { return ((inter_case<0) ? -inter_case : 0); }

private:
  G4InuclParticle* bullet;
  G4InuclParticle* target;

  G4int inter_case;
};

#endif // G4INTERACTION_CASE_HH 


