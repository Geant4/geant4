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
// $Id: G4InteractionCase.cc 66241 2012-12-13 18:34:42Z gunter $
//
// 20100518  M. Kelsey -- Move code from Colliders' "bulletTargetSetter()"
//		to setBulletTarget().

#include "G4InteractionCase.hh"
#include "G4InuclParticle.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"


// Evaluate concrete types of input particles and assign bullet, target

void G4InteractionCase::set(G4InuclParticle* part1, 
			    G4InuclParticle* part2) {
  clear();		// Reset everything in case of failure

  // See which one of the two (or both) is a nucleus
  G4InuclNuclei* nucl1 = dynamic_cast<G4InuclNuclei*>(part1);
  G4InuclNuclei* nucl2 = dynamic_cast<G4InuclNuclei*>(part2);

  G4InuclElementaryParticle* had1 = dynamic_cast<G4InuclElementaryParticle*>(part1);
  G4InuclElementaryParticle* had2 = dynamic_cast<G4InuclElementaryParticle*>(part2);

  if (nucl1 && nucl2) { 	// Nuclear collision, lighter is projectile
    inter_case = -2;
    if (nucl2->getA() >= nucl1->getA()) {
      bullet = part1;
      target = part2;
    } else {
      bullet = part2;
      target = part1;
    } 
  } else if (nucl1 || nucl2) {	// Hadron on nucleus, hadron projectile
    inter_case = -1;
    if (nucl1 && had2) {
      bullet = part2;
      target = part1;
    } else {
      bullet = part1;
      target = part2;
    }
  } else if (had1 && had2) {	// Hadron-hadron interaction, order irrelevant
    inter_case = had1->type() * had2->type();
    bullet = part1;
    target = part2;
  }
}
