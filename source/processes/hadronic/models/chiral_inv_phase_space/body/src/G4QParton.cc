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
// $Id: G4QParton.cc,v 1.6 2009-07-02 07:17:09 mkossov Exp $
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

G4QParton::G4QParton()
{
  // CHIPS is working only with u, d, and s quarks (SU(3)xSU(3)) (no gluons! M.K.)
  // Random Flavor/Colour/Spin definition for default constructor (with .3 s-suppresion)
  PDGencoding=(G4int)(2.3*G4UniformRand())+1; //@@ What about antiquarks? (M.K.)
  theDefinition=G4ParticleTable::GetParticleTable()->FindParticle(PDGencoding);
#ifdef debug
  G4cout<<"G4QParton::DefConstructer: Type = "<<theDefinition->GetParticleSubType()<<G4endl;
#endif
  // random colour (1,2,3)=(R,G,B) for quarks and (-1,-2,-3)=(aR,aG,aB) for anti-quarks
  theColour = (G4int)(3*G4UniformRand())+1;
  if(theColour>3) theColour = 3;
  theIsoSpinZ = theDefinition->GetPDGIsospin3();
  theSpinZ = (G4int)(2*G4UniformRand()) - 0.5;
}
G4QParton::G4QParton(G4int PDGcode)
{
  PDGencoding=PDGcode;
  theX = 0;
  theDefinition=G4ParticleTable::GetParticleTable()->FindParticle(PDGencoding);
#ifdef debug
  G4cout<<"G4QParton::PDGConstructer: Type = "<<theDefinition->GetParticleSubType()<<G4endl;
#endif
  G4int aPDG=std::abs(PDGcode);
  //___________quarks__ 2 possible codes for gluons _ Condition for di-quarks
  if(!aPDG || (aPDG>3 && PDGcode!=9 && PDGcode!=21 &&
          (aPDG>3303||aPDG<1103||aPDG%100>3)) || theDefinition==0)
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
  // ColourByRandom (-11,-12,-13,-21,...,-33)=(RRbar,RGbar,RBbar,...,BBbar) for gluons
  else if(name == "gluons") theColour = -(RGB*10 + (G4int)(3*G4UniformRand())+1);
  else
  {
    G4cerr<<"***G4QParton::Constructor: not quark/diquark/gluon = "
          <<theDefinition->GetParticleType()<<G4endl;
    G4Exception("G4QParton::Constructor:","72",FatalException,"WrongParton");
  }
  // isospin-z from PDG-encoded isospin-z for 
  // quarks, anti-quarks, di-quarks, and anti-di-quarks:
  if(name == "quarks" || name == "diquarks") theIsoSpinZ = theDefinition->GetPDGIsospin3();
       // isospin-z choosen at random from PDG-encoded isospin for gluons (should be zero):
  else
  {
    G4int thisPDGiIsospin=theDefinition->GetPDGiIsospin();
    if (thisPDGiIsospin == 0) theIsoSpinZ = 0;
         //@@ ? M.K.
    else theIsoSpinZ=((G4int)((thisPDGiIsospin+1)*G4UniformRand())) - thisPDGiIsospin*0.5;
  }
  //
  // spin-z choosen at random from PDG-encoded spin:
  //
  G4int thisPDGiSpin=theDefinition->GetPDGiSpin();
  if(thisPDGiSpin == 0) theSpinZ = 0;
  else theSpinZ = (G4int)((thisPDGiSpin+1)*G4UniformRand())-thisPDGiSpin*0.5;;
}

G4QParton::G4QParton(const G4QParton &right)
{
  PDGencoding = right.PDGencoding;
  theMomentum = right.theMomentum;
  thePosition = right.thePosition;
  theX = right.theX;
  theDefinition = right.theDefinition;
#ifdef debug
  G4cout<<"G4QParton::RCopyConstructer: Type="<<theDefinition->GetParticleSubType()<<G4endl;
#endif
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
#ifdef debug
  G4cout<<"G4QParton::PCopyConstructer: Type="<<theDefinition->GetParticleSubType()<<G4endl;
#endif
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
#ifdef debug
    G4cout<<"G4QParton::=Constructer: Type = "<<theDefinition->GetParticleSubType()<<G4endl;
#endif
    theColour = right.theColour;
    theIsoSpinZ = right.theIsoSpinZ;
    theSpinZ = right.theSpinZ; 
  }
  return *this;
}

G4QParton::~G4QParton()
{
  //G4cout << "G4QParton::~G4QParton(): this = "<<this <<G4endl;
  //G4cout << "break here"<<this <<G4endl;
}
// QGS x+/x- logic of the Energy and Pz calculation
void G4QParton::DefineMomentumInZ(G4double aLightConeMomentum, G4bool aDirection)
{
  G4double Mass = GetMass();                       // Should be zero for u,d
  G4LorentzVector a4Momentum = Get4Momentum();
  aLightConeMomentum*=theX;
  G4double TransverseMass2 = sqr(a4Momentum.px()) + sqr(a4Momentum.py()) + sqr(Mass);
  a4Momentum.setPz(0.5*(aLightConeMomentum - TransverseMass2/aLightConeMomentum) *
                                                                      (aDirection? 1: -1));
  a4Momentum.setE( 0.5*(aLightConeMomentum + TransverseMass2/aLightConeMomentum));
  Set4Momentum(a4Momentum);
}  
