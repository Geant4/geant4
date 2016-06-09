//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4Nucleon.cc,v 1.6 2002/12/12 19:17:58 gunter Exp $
// GEANT4 tag $Name: geant4-05-01 $
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
: theBindingE(0.) , theParticleType(NULL), theSplitableHadron(NULL)
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

#include "g4std/iostream"
G4std::ostream & operator << (G4std::ostream &s, const G4Nucleon& nucleon)
{
//	s<< nucleon.GetDefinition()->GetParticleName() 
//	 << "  is " << nucleon.AreYouHit() ? " " : "not" 
//	 << " hit. Momentum/position:" << G4endl;
	s<< "  momentum : " << nucleon.Get4Momentum() << G4endl;
	s<< "  position : " << nucleon.GetPosition() ;
	return s;
}	  
