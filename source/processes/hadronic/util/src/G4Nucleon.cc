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
//
//
#include "G4Nucleon.hh"

// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4Nucleon ----------------
//             by Gunter Folger, May 1998.
//       class for a nucleon (inside a 3D Nucleus)
// ------------------------------------------------------------

G4Nucleon::G4Nucleon()
: theBindingE(0.) , theParticleType(0), theSplitableHadron(0)
{}

G4Nucleon::~G4Nucleon()
{
}

void G4Nucleon::Boost(const G4LorentzVector & aMomentum)
{
//   see e.g. CERNLIB short writeup U101 for the algorithm
	G4double mass=aMomentum.mag();
	G4double factor=
	    ( theMomentum.vect()*aMomentum.vect()/(aMomentum.e()+mass) - theMomentum.e() ) / mass;

	theMomentum.setE(1/mass*theMomentum.dot(aMomentum));
	theMomentum.setVect(factor*aMomentum.vect() + theMomentum.vect());
}

#include <iostream>
std::ostream & operator << (std::ostream &stream, const G4Nucleon& nucleon)
{
//	stream<< nucleon.GetDefinition()->GetParticleName()
//	 << "  is " << nucleon.AreYouHit() ? " " : "not" 
//	 << " hit. Momentum/position:" << G4endl;
	stream<< "  momentum : " << nucleon.Get4Momentum() << G4endl;
	stream<< "  position : " << nucleon.GetPosition() ;
	return stream;
}	  
