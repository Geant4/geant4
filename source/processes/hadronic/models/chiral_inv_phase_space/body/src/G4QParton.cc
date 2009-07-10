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
// $Id: G4QParton.cc,v 1.7 2009-07-10 16:42:57 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4QParton ----------------
//             by Mikhail Kossov, Oct 2006.
//   class for Parton (inside a string) used by Parton String Models
//   For comparison mirror member functions are taken from G4 class:
//   G4Parton
// ------------------------------------------------------------------------
// Short description: The Quark-Gluon String consists of the partons, which
// are quarks and some times gluons.
// ------------------------------------------------------------------------

//#define debug

#include "G4QParton.hh"

G4QParton::G4QParton() // By default creates only quarks (not di-quarks)
{
  // CHIPS is working only with u, d, and s quarks (SU(3)xSU(3)) (no gluons! M.K.)
  // Random Flavor/Colour/Spin definition for default constructor (with .3 s-suppresion)
  PDGencoding=(G4int)(2.3*G4UniformRand())+1; //@@ Additional parameter of s/u (M.K.)
  theType=1;
#ifdef debug
  G4cout<<"....G4QParton::DefConstructer: PDG = "<<PDGencoding<<", Type="<<theType<<G4endl;
#endif
  // random colour (1,2,3)=(R,G,B) for quarks and (-1,-2,-3)=(aR,aG,aB) for anti-quarks
  theColour = (G4int)(3*G4UniformRand())+1;
  if(theColour>3) theColour = 3;                   // Should never happend
  theSpinZ = (G4int)(2*G4UniformRand()) - 0.5;
  // Default definition (initialization)
  theX = 0.;
  thePosition=G4ThreeVector(0.,0.,0.);
  theMomentum=G4LorentzVector(0.,0.,0.,0.);
}

G4QParton::G4QParton(G4int PDGcode)
{
  SetPDGCode(PDGcode);
  // Default definition (initialization)
  theX = 0.;
  thePosition=G4ThreeVector(0.,0.,0.);
  theMomentum=G4LorentzVector(0.,0.,0.,0.);
}

G4QParton::G4QParton(const G4QParton &right)
{
  PDGencoding = right.PDGencoding;
  theType = right.theType;
  theMomentum = right.theMomentum;
  thePosition = right.thePosition;
  theX = right.theX;
  theColour = right.theColour;
  theSpinZ = right.theSpinZ;
#ifdef debug
  G4cout<<"G4QParton::RCopyConstructer: PDG="<<PDGencoding<<", Col="<<theColour<<", Sz="
        <<theSpinZ<<G4endl;
#endif
}

G4QParton::G4QParton(const G4QParton* right)
{
  PDGencoding   = right->PDGencoding;
  theType       = right->theType;
  theMomentum   = right->theMomentum;
  thePosition   = right->thePosition;
  theX          = right->theX;
  theColour     = right->theColour;
  theSpinZ      = right->theSpinZ;
#ifdef debug
  G4cout<<"G4QParton::PCopyConstructer: PDG="<<PDGencoding<<", Col="<<theColour<<", Sz="
        <<theSpinZ<<G4endl;
#endif
}

const G4QParton& G4QParton::operator=(const G4QParton &right)
{
  if(this != &right)                          // Beware of self assignment
  {
    PDGencoding=right.GetPDGCode();
    theType=right.GetType();
    theMomentum=right.Get4Momentum();
    thePosition=right.GetPosition();
    theX = right.theX;
    theColour = right.theColour;
    theSpinZ = right.theSpinZ; 
#ifdef debug
    G4cout<<"G4QParton::=Constructer: PDG="<<PDGencoding<<", Col="<<theColour<<", Sz="
          <<theSpinZ<<G4endl;
#endif
  }
  return *this;
}

G4QParton::~G4QParton() {}

// Redefine the parton nature without changing x, 4Mom, Pos etc.
void G4QParton::SetPDGCode(G4int PDGcode)
{
  PDGencoding=PDGcode;
  G4int aPDG=std::abs(PDGcode);
  if(aPDG < 3304 && aPDG > 1100 && aPDG%100 < 4) theType=2; // di-quark
  else if(aPDG && aPDG<4)                        theType=1; // quark
  else if(aPDG==9 || aPDG==21)                   theType=0; // gluon
  else
  {
    G4cerr<<"***G4QParton::SetPDGCode: wrong gluon/quark/diquark PDG="<<PDGcode<<G4endl;
    G4Exception("G4QParton::SetPDGCode:","72",FatalException,"WrongPartonPDG");
  }
#ifdef debug
  G4cout<<"....G4QParton::SetPDGCode: PDG = "<<PDGcode<<", Type="<<theType<<G4endl;
#endif
  //
  // colour by random in (1,2,3)=(R,G,B) for quarks and 
  //                  in (-1,-2,-3)=(Rbar,Gbar,Bbar) for anti-quarks:
  G4int RGB=(G4int)(3*G4UniformRand())+1;
  if(theType==1)
  {
    if(PDGcode>0) theColour = RGB;
    else          theColour =-RGB;
  }
  // colour by random in (-1,-2,-3)=(Rbar,Gbar,Bbar)=(GB,RB,RG) for di-quarks and
  //                  in (1,2,3)=(R,G,B)=(GB,RB,RG) for anti-di-quarks:
  else if(theType==2)
  {
    if(PDGcode>0) theColour =-RGB;
    else          theColour = RGB;
  }
  // ColourByRandom (-11,-12,-13,-21,...,-33)=(RRbar,RGbar,RBbar,...,BBbar) for gluons
  else theColour = -(RGB*10 + (G4int)(3*G4UniformRand())+1);
  //
  // spin-z choosen at random from PDG-encoded spin:
  //
  G4double            dPDGSpin=1.;        // Quark 2S
  if     (theType==0) dPDGSpin=2.;        // Gluon 2S
  else if(theType==2) dPDGSpin=aPDG%10-1; // Di-quark 2S
  theSpinZ = (G4int)((dPDGSpin+1)*G4UniformRand())-dPDGSpin/2;
}

// QGS x+/x- logic of the Energy and Pz calculation
void G4QParton::DefineMomentumInZ(G4double aLightConeMomentum, G4bool aDirection)
{
  G4LorentzVector a4Momentum = Get4Momentum();
  aLightConeMomentum*=theX;
  G4double TransverseMass2 = sqr(a4Momentum.px()) + sqr(a4Momentum.py());
  a4Momentum.setPz(0.5*(aLightConeMomentum - TransverseMass2/aLightConeMomentum) *
                                                                      (aDirection? 1: -1));
  a4Momentum.setE( 0.5*(aLightConeMomentum + TransverseMass2/aLightConeMomentum));
  Set4Momentum(a4Momentum);
}  
