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
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4Parton ----------------
//             by Gunter Folger, June 1998.
//       class for Parton (inside a string) used by Parton String Models
// ------------------------------------------------------------

#include "G4Parton.hh"
#include "G4HadronicException.hh"

G4Parton::G4Parton(G4int PDGcode)
{
	PDGencoding=PDGcode;
	theX = 0;
	theDefinition=G4ParticleTable::GetParticleTable()->FindParticle(PDGencoding);
	if (theDefinition == NULL)
	{
	  G4cout << "Encoding = "<<PDGencoding<<G4endl;
	  G4String text = "G4Parton::GetDefinition(): Encoding not in particle table";
	  throw G4HadronicException(__FILE__, __LINE__, text);
	}
	//
	// colour by random in (1,2,3)=(R,G,B) for quarks and 
	//                  in (-1,-2,-3)=(Rbar,Gbar,Bbar) for anti-quarks:
  //
	if (theDefinition->GetParticleType() == "quarks") {
		theColour = ((G4int)(3.*G4UniformRand())+1)*(std::abs(PDGencoding)/PDGencoding) ;
	}
	//
	// colour by random in (-1,-2,-3)=(Rbar,Gbar,Bbar)=(GB,RB,RG) for di-quarks and
	//                  in (1,2,3)=(R,G,B)=(GB,RB,RG) for anti-di-quarks:
  //
	else if (theDefinition->GetParticleType() == "diquarks") {
		theColour = -((G4int)(3.*G4UniformRand())+1)*(std::abs(PDGencoding)/PDGencoding);
	}
	//
	// colour by random in (-11,-12,...,-33)=(RRbar,RGbar,RBbar,...,BBbar) for gluons:
  //
	else if (theDefinition->GetParticleType() == "gluons") {
		theColour = -(((G4int)(3.*G4UniformRand())+1)*10 + ((G4int)(3.*G4UniformRand())+1));
	}
	else {
	  G4cout << "Encoding = "<<PDGencoding<<G4endl;
	  G4String text = "G4Parton::GetDefinition(): Particle is not a parton";
	  throw G4HadronicException(__FILE__, __LINE__, text);
	}
	//  
	// isospin-z from PDG-encoded isospin-z for 
	// quarks, anti-quarks, di-quarks, and anti-di-quarks:
  //
	if ((theDefinition->GetParticleType() == "quarks") || (theDefinition->GetParticleType() == "diquarks")){
		theIsoSpinZ = theDefinition->GetPDGIsospin3();
	}
	//
  // isospin-z choosen at random from PDG-encoded isospin for gluons (should be zero):
	//
	else {
		G4int thisPDGiIsospin=theDefinition->GetPDGiIsospin();
		if (thisPDGiIsospin == 0) {
			theIsoSpinZ = 0;
		}
		else {
			theIsoSpinZ = ((G4int)((thisPDGiIsospin+1)*G4UniformRand()))-thisPDGiIsospin*0.5;
		}
	}
	//
	// spin-z choosen at random from PDG-encoded spin:
	//
	G4int thisPDGiSpin=theDefinition->GetPDGiSpin();
	if (thisPDGiSpin == 0) {
		theSpinZ = 0;
	}
	else {
		G4int rand=((G4int)((thisPDGiSpin+1)*G4UniformRand()));
		theSpinZ = rand-thisPDGiSpin*0.5;;
	}
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

G4Parton & G4Parton::operator=(const G4Parton &right)
{
   if (this != &right)
   {
      PDGencoding=right.GetPDGcode();
      theMomentum=right.Get4Momentum();
      thePosition=right.GetPosition();
      theX = right.theX;
      theDefinition = right.theDefinition;
      theColour = right.theColour;
      theIsoSpinZ = right.theIsoSpinZ;
      theSpinZ = right.theSpinZ;
   }

	return *this;
}

G4Parton::~G4Parton()
{
//  cout << "G4Parton::~G4Parton(): this = "<<this <<endl;
//  cout << "break here"<<this <<endl;
}

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

void G4Parton::DefineMomentumInZ(G4double aLightConeMomentum,G4double aLightConeE, G4bool aDirection)
{
	G4double Mass = GetMass();
	G4LorentzVector a4Momentum = Get4Momentum();
	aLightConeMomentum*=theX;
	aLightConeE*=theX;
	G4double TransverseMass2 = sqr(a4Momentum.px()) + sqr(a4Momentum.py()) + sqr(Mass);
	a4Momentum.setPz(0.5*(aLightConeMomentum - aLightConeE - TransverseMass2/aLightConeMomentum)*(aDirection? 1: -1)); 
	a4Momentum.setE( 0.5*(aLightConeMomentum + aLightConeE + TransverseMass2/aLightConeMomentum));
	Set4Momentum(a4Momentum);
}  
