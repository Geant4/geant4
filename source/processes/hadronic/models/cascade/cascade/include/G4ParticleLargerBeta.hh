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
// $Id: G4ParticleLargerBeta.hh 70080 2013-05-22 22:46:28Z mkelsey $
//
// 20100112  M. Kelsey -- Add additional operator() which uses pointers,
//		also fix bug which returns wrong result
// 20110415  M. Kelsey -- Add additional operator() for G4CascadParticle
// 20130522  M. Kelsey -- Long-standing error: should be ">", not ">="

#ifndef G4ParticleLargerBeta_h
#define G4ParticleLargerBeta_h

#include "G4CascadParticle.hh"
#include "G4InuclElementaryParticle.hh"

#ifdef G4CASCADE_DEBUG_SORT
#include "G4ios.hh"
#endif

class G4ParticleLargerBeta {
public:
  G4bool operator() (const G4InuclElementaryParticle& part1,
		     const G4InuclElementaryParticle& part2) {
#ifdef G4CASCADE_DEBUG_SORT
    G4cout << "part1 @ " << &part1 << ": " << part1
	   << "part2 @ " << &part2 << ": " << part2
	   << G4endl;
#endif
    return (part1.getMomModule()/part1.getEnergy() >
	    part2.getMomModule()/part2.getEnergy()
	    );
  }
 
  G4bool operator() (const G4InuclElementaryParticle* part1,
		     const G4InuclElementaryParticle* part2) {
    return (part1 && part2 && operator()(*part1, *part2));
  }

  G4bool operator() (const G4CascadParticle& part1,
		     const G4CascadParticle& part2) {
    return (operator()(part1.getParticle(), part2.getParticle()));
  }

  G4bool operator() (const G4CascadParticle* part1,
		     const G4CascadParticle* part2) {
    return (part1 && part2 && operator()(*part1, *part2));
  }
};

#endif	/* G4ParticleLargerBeta_h */
