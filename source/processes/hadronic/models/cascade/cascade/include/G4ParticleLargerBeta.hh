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
// $Id: G4ParticleLargerBeta.hh,v 1.7 2010-06-25 09:43:44 gunter Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100112  M. Kelsey -- Add additional operator() which uses pointers,
//		also fix bug which returns wrong result

#ifndef G4ParticleLargerBeta_h
#define G4ParticleLargerBeta_h

#include "G4InuclElementaryParticle.hh"

class G4ParticleLargerBeta {
public:
  G4bool operator() (const G4InuclElementaryParticle& part1,
		     const G4InuclElementaryParticle& part2) {
    return (part1.getMomModule()/part1.getEnergy() >=
	    part2.getMomModule()/part2.getEnergy()
	    );
  }
 
  G4bool operator() (const G4InuclElementaryParticle* part1,
		     const G4InuclElementaryParticle* part2) {
    return (part1 && part2 && operator()(*part1, *part2));
  }
};

#endif	/* G4ParticleLargerBeta_h */
