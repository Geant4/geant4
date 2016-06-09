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
// $Id: G4QParton.cc,v 1.3 2006/12/12 11:02:22 mkossov Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4QParton ----------------
//             by Mikhail Kossov, Oct 2006.
//   class for Parton (inside a string) used by Parton String Models
//   For comparison mirror member functions are taken from G4 class:
//   G4Parton
// ------------------------------------------------------------

#include "G4QParton.hh"

G4QParton::G4QParton(G4int PDGcode)
{
	 PDGencoding=PDGcode;
	 theX = 0;
	 theDefinition=G4ParticleTable::GetParticleTable()->FindParticle(PDGencoding);
  G4int aPDG=std::abs(PDGcode);
	 //___________quarks__ 2 possible codes for gluons _ Condition for di-quarks
	 if(!aPDG || (aPDG>3 && PDGcode!=9 && PDGcode!=21 && (aPDG>3303||aPDG<1103||aPDG%100>3))
     || theDefinition==0)
	 {
	 		G4cerr<<"***G4QParton::Constructor: wrong quark/diquark PDG="<<PDGcode<<G4endl;
    G4Exception("G4QParton::Constructor:","72",FatalException,"WrongPartonPDG");
	 }
	 //
	 // colour by random in (1,2,3)=(R,G,B) for quarks and 
	 //                  in (-1,-2,-3)=(Rbar,Gbar,Bbar) for anti-quarks:
  G4int RGB=(G4int)(3*G4UniformRand())+1;
	 G4String name=theDefinition->GetParticleType();
	 if(name == "quarks")
	 {
	  	if(PDGcode>0) theColour = RGB;
    else          theColour =-RGB;
  }
	 // colour by random in (-1,-2,-3)=(Rbar,Gbar,Bbar)=(GB,RB,RG) for di-quarks and
	 //                  in (1,2,3)=(R,G,B)=(GB,RB,RG) for anti-di-quarks:
	 else if(name == "diquarks")
  {
	  	if(PDGcode>0) theColour =-RGB;
    else          theColour = RGB;
	 }
	 // colour by random in (-11,-12,-13,-21,...,-33)=(RRbar,RGbar,RBbar,...,BBbar) for gluons
	 else if(name == "gluons") theColour = -(RGB*10 + (G4int)(3*G4UniformRand())+1);
	 else
  {
	 		G4cerr<<"***G4QParton::Constructor: not quark/diquark/gluon = "
          <<theDefinition->GetParticleType()<<G4endl;
    G4Exception("G4QParton::Constructor:","72",FatalException,"WrongParton");
	 }
	 // isospin-z from PDG-encoded isospin-z for 
	 // quarks, anti-quarks, di-quarks, and anti-di-quarks:
	 if (name == "quarks" || name == "diquarks")theIsoSpinZ = theDefinition->GetPDGIsospin3();
  // isospin-z choosen at random from PDG-encoded isospin for gluons (should be zero):
	 else
  {
		  G4int thisPDGiIsospin=theDefinition->GetPDGiIsospin();
		  if (thisPDGiIsospin == 0) theIsoSpinZ = 0;
    //@@ ? M.K.
		  else	theIsoSpinZ = ((G4int)((thisPDGiIsospin+1)*G4UniformRand()))-thisPDGiIsospin*0.5;
	 }
	 //
	 // spin-z choosen at random from PDG-encoded spin:
	 //
	 G4int thisPDGiSpin=theDefinition->GetPDGiSpin();
	 if(thisPDGiSpin == 0)	theSpinZ = 0;
	 else	theSpinZ = (G4int)((thisPDGiSpin+1)*G4UniformRand())-thisPDGiSpin*0.5;;
}

G4QParton::G4QParton(const G4QParton &right)
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

G4QParton::G4QParton(const G4QParton* right)
{
  PDGencoding   = right->PDGencoding;
  theMomentum   = right->theMomentum;
  thePosition   = right->thePosition;
  theX          = right->theX;
  theDefinition = right->theDefinition;
  theColour     = right->theColour;
  theIsoSpinZ   = right->theIsoSpinZ;
  theSpinZ      = right->theSpinZ;
}

const G4QParton& G4QParton::operator=(const G4QParton &right)
{
  if(this != &right)                          // Beware of self assignment
  {
    PDGencoding=right.GetPDGCode();
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

G4QParton::~G4QParton()
{
//  cout << "G4QParton::~G4QParton(): this = "<<this <<endl;
//  cout << "break here"<<this <<endl;
}

void G4QParton::DefineMomentumInZ(G4double aLightConeMomentum, G4bool aDirection)
{
	G4double Mass = GetMass();
	G4LorentzVector a4Momentum = Get4Momentum();
	aLightConeMomentum*=theX;
	G4double TransverseMass2 = sqr(a4Momentum.px()) + sqr(a4Momentum.py()) + sqr(Mass);
	a4Momentum.setPz(0.5*(aLightConeMomentum - TransverseMass2/aLightConeMomentum)*(aDirection? 1: -1)); 
	a4Momentum.setE( 0.5*(aLightConeMomentum + TransverseMass2/aLightConeMomentum));
	Set4Momentum(a4Momentum);
}  
