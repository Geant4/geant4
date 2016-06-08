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
// $Id: G4Parton.cc,v 1.7.8.2 2001/06/28 20:20:03 gunter Exp $
// GEANT4 tag $Name:  $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4Parton ----------------
//             by Gunter Folger, June 1998.
//       class for Parton (inside a string) used by Parton String Models
// ------------------------------------------------------------

#include "G4Parton.hh"

G4Parton::G4Parton(G4int PDGcode)
{
	PDGencoding=PDGcode;
	theX = 0;
	theDefinition=G4ParticleTable::GetParticleTable()->FindParticle(PDGencoding);;
	if (theDefinition == NULL)
	{
	  G4cout << "Encoding = "<<PDGencoding<<G4endl;
	  G4Exception("G4Parton::GetDefinition(): Encoding not in particle table");
	}
	theColour = 1;
	theIsoSpinZ = 0.5;
	theSpinZ = 0.5;
}

G4Parton::G4Parton(const G4Parton &right)
{
	PDGencoding = right.PDGencoding;
	theMomentum = right.theMomentum;
	thePosition = right.thePosition;
	theX = right.theX;
	theDefinition = right.theDefinition;
	theColour = right.theColour;
	theIsoSpinZ = right.theIsoSpinZ;
	theSpinZ = right.theSpinZ;
}

const G4Parton & G4Parton::operator=(const G4Parton &right)
{
	PDGencoding=right.GetPDGcode();
	theMomentum=right.Get4Momentum();
	thePosition=right.GetPosition();
	theX = right.theX;
	theDefinition = right.theDefinition;
	theColour = right.theColour;
	theIsoSpinZ = right.theIsoSpinZ;
	theSpinZ = right.theSpinZ;
		
	return *this;
}

G4Parton::~G4Parton()
{}

void G4Parton::DefineMomentumInZ(G4double aLightConeMomentum, G4bool aDirection)
{
	G4double Mass = GetMass();
	G4LorentzVector a4Momentum = Get4Momentum();
	aLightConeMomentum*=theX;
	G4double TransverseMass2 = sqr(a4Momentum.px()) + sqr(a4Momentum.py()) + sqr(Mass);
	a4Momentum.setPz(0.5*(aLightConeMomentum - TransverseMass2/aLightConeMomentum)*(aDirection? 1: -1)); 
	a4Momentum.setE( 0.5*(aLightConeMomentum + TransverseMass2/aLightConeMomentum));
	Set4Momentum(a4Momentum);
}  
