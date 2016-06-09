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
// $Id$
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
  PGGCode=(G4int)(2.3*G4UniformRand())+1; //@@ Additional parameter of s/u (M.K.)
  theType=1;
#ifdef debug
  G4cout<<"....G4QParton::DefConstructer: PDG = "<<PGGCode<<", Type="<<theType<<G4endl;
#endif
  // random colour (1,2,3)=(R,G,B) for quarks and (-1,-2,-3)=(aR,aG,aB) for anti-quarks
  theColour = (G4int)(3*G4UniformRand())+1;
  if(theColour>3) theColour = 3;                   // Should never happend
  theSpinZ = (G4int)(2*G4UniformRand()) - 0.5;
  QCont = G4QContent(0,0,0,0,0,0);  
  // Default definition (initialization)
  theX = 0.;
  thePosition=G4ThreeVector(0.,0.,0.);
  theMomentum=G4LorentzVector(0.,0.,0.,0.);
}

G4QParton::G4QParton(G4int PDG)
{
  SetPDGCode(PDG);
  // Default definition (initialization)
  theX = 0.;
  thePosition=G4ThreeVector(0.,0.,0.);
  theMomentum=G4LorentzVector(0.,0.,0.,0.);
}

G4QParton::G4QParton(const G4QParton &right)
{
  PGGCode     = right.PGGCode;
  QCont       = right.QCont;
  theType     = right.theType;
  theMomentum = right.theMomentum;
  thePosition = right.thePosition;
  theX        = right.theX;
  theColour   = right.theColour;
  theSpinZ    = right.theSpinZ;
#ifdef debug
  G4cout<<"G4QParton::RCopyConstructer: PDG="<<PGGCode<<", Col="<<theColour<<", Sz="
        <<theSpinZ<<G4endl;
#endif
}

G4QParton::G4QParton(const G4QParton* right)
{
  PGGCode       = right->PGGCode;
  QCont         = right->QCont;
  theType       = right->theType;
  theMomentum   = right->theMomentum;
  thePosition   = right->thePosition;
  theX          = right->theX;
  theColour     = right->theColour;
  theSpinZ      = right->theSpinZ;
#ifdef debug
  G4cout<<"G4QParton::PCopyConstructer: PDG="<<PGGCode<<", Col="<<theColour<<", Sz="
        <<theSpinZ<<G4endl;
#endif
}

const G4QParton& G4QParton::operator=(const G4QParton &right)
{
  if(this != &right)                          // Beware of self assignment
  {
    PGGCode     = right.GetPDGCode();
    QCont       = right.QCont;
    theType     = right.GetType();
    theMomentum = right.Get4Momentum();
    thePosition = right.GetPosition();
    theX        = right.theX;
    theColour   = right.theColour;
    theSpinZ    = right.theSpinZ; 
#ifdef debug
    G4cout<<"G4QParton::=Constructer: PDG="<<PGGCode<<", Col="<<theColour<<", Sz="
          <<theSpinZ<<G4endl;
#endif
  }
  return *this;
}

G4QParton::~G4QParton() {}

// Redefine the parton nature without changing x, 4Mom, Pos etc.
void G4QParton::SetPDGCode(G4int PDG)
{
  PGGCode=PDG;
  G4int aPDG=std::abs(PDG);
  if(aPDG < 3304 && aPDG > 1100 && aPDG%100 < 4) // di-quark
  {
    theType=2;
    G4int cPDG=aPDG/100;
    if(PDG>0)
    {
      if     (cPDG==11) QCont=G4QContent(2,0,0,0,0,0);   // dd
      else if(cPDG==21) QCont=G4QContent(1,1,0,0,0,0);   // ud
      else if(cPDG==22) QCont=G4QContent(0,2,0,0,0,0);   // uu
      else if(cPDG==31) QCont=G4QContent(1,0,1,0,0,0);   // sd
      else if(cPDG==32) QCont=G4QContent(0,1,1,0,0,0);   // su
      else if(cPDG==33) QCont=G4QContent(0,0,2,0,0,0);   // ss
      else
      {
        G4cerr<<"***G4QParton::SetPDGCode: bad di-quark PDG="<<PDG<<G4endl;
        G4Exception("G4QParton::SetPDGCode:","72",FatalException,"Not SU(3) DiQuark");
      }
    }
    else
    {
      if     (cPDG==11) QCont=G4QContent(0,0,0,2,0,0);   // anti-dd
      else if(cPDG==21) QCont=G4QContent(0,0,0,1,1,0);   // anti-ud
      else if(cPDG==22) QCont=G4QContent(0,0,0,0,2,0);   // anti-uu
      else if(cPDG==31) QCont=G4QContent(0,0,0,1,0,1);   // anti-sd
      else if(cPDG==32) QCont=G4QContent(0,0,0,0,1,1);   // anti-su
      else if(cPDG==33) QCont=G4QContent(0,0,0,0,0,2);   // anti-ss
      else
      {
        G4cerr<<"***G4QParton::SetPDGCode: bad anti-di-quark PDG="<<PDG<<G4endl;
        G4Exception("G4QParton::SetPDGCode:","72",FatalException,"Not SU(3) AntiDiQuark");
      }
    }
  }
  else if(aPDG && aPDG<4)                        // quark
  {
    theType=1;
    if(PDG>0)
    {
      if     (PDG==1) QCont=G4QContent(1,0,0,0,0,0);   // d
      else if(PDG==2) QCont=G4QContent(0,1,0,0,0,0);   // u
      else if(PDG==3) QCont=G4QContent(0,0,1,0,0,0);   // s
      else
      {
        G4cerr<<"***G4QParton::SetPDGCode: bad quark PDG="<<PDG<<G4endl;
        G4Exception("G4QParton::SetPDGCode:","72",FatalException,"Not SU(3) Quark");
      }
    }
    else
    {
      if     (PDG==-1) QCont=G4QContent(0,0,0,1,0,0);  // anti-d
      else if(PDG==-2) QCont=G4QContent(0,0,0,0,1,0);  // anti-u
      else if(PDG==-3) QCont=G4QContent(0,0,0,0,0,1);  // anti-s
      else
      {
        G4cerr<<"***G4QParton::SetPDGCode: bad anti-quark PDG="<<PDG<<G4endl;
        G4Exception("G4QParton::SetPDGCode:","72",FatalException,"Not SU(3) Anti-Quark");
      }
    }
  }
  else if(aPDG==9 || aPDG==21)                   // gluon
  {
    theType=0;
    QCont=G4QContent(0,0,0,0,0,0);
  }
  else
  {
    G4cerr<<"***G4QParton::SetPDGCode: wrong gluon/quark/diquark PDG="<<PDG<<G4endl;
    G4Exception("G4QParton::SetPDGCode:","72",FatalException,"WrongPartonPDG");
  }
#ifdef debug
  G4cout<<"....G4QParton::SetPDGCode: PDG = "<<PDG<<", Type="<<theType<<G4endl;
#endif
  //
  // colour by random in (1,2,3)=(R,G,B) for quarks and 
  //                  in (-1,-2,-3)=(Rbar,Gbar,Bbar) for anti-quarks:
  G4int RGB=(G4int)(3*G4UniformRand())+1;
  if(theType==1)
  {
    if(PDG>0) theColour = RGB;
    else      theColour =-RGB;
  }
  // colour by random in (-1,-2,-3)=(Rbar,Gbar,Bbar)=(GB,RB,RG) for di-quarks and
  //                  in (1,2,3)=(R,G,B)=(GB,RB,RG) for anti-di-quarks:
  else if(theType==2)
  {
    if(PDG>0) theColour =-RGB;
    else      theColour = RGB;
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

// Reduce DiQ-aDiQ to Q-aQ (true if succeeded). General function of the QPartons operations
G4bool G4QParton::ReduceDiQADiQ(G4QParton* d1, G4QParton* d2)
{
  G4bool result=false;
  G4int sPDG=d1->GetPDGCode();
  G4int nPDG=d2->GetPDGCode();
#ifdef debug
  G4cout<<"G4QParton::ReduceDiQADiQ: **Called** LPDG="<<sPDG<<", RPDG="<<nPDG<<G4endl;
#endif
  G4int        qPDG=sPDG;
  if(qPDG<-99) qPDG=(-qPDG)/100;
  else         qPDG/=100;
  G4int        dPDG=nPDG;
  if(dPDG<-99) dPDG=(-dPDG)/100;
  else         dPDG/=100;
  G4int L1=qPDG/10;
  G4int L2=qPDG%10;
  G4int R1=dPDG/10;
  G4int R2=dPDG%10;
  if(L1==R1 || L1==R2 || L2==R1 || L2==R2) // Annihilation condition
  {
    if     (L1==R1)
    {
      if(sPDG>0) sPDG=L2;
      else       sPDG=-L2;
      if(nPDG>0) nPDG=R2;
      else       nPDG=-R2;
#ifdef debug
      G4cout<<"G4QParton::ReDiQADiQ:L2="<<L2<<",R2="<<R2<<",L="<<sPDG<<",R="<<nPDG<<G4endl;
#endif
    }
    else if(L1==R2)
    {
      if(sPDG>0) sPDG=L2;
      else       sPDG=-L2;
      if(nPDG>0) nPDG=R1;
      else       nPDG=-R1;
#ifdef debug
      G4cout<<"G4QParton::ReDiQADiQ:L2="<<L2<<",R1="<<R1<<",L="<<sPDG<<",R="<<nPDG<<G4endl;
#endif
    }
    else if(L2==R1)
    {
      if(sPDG>0) sPDG=L1;
      else       sPDG=-L1;
      if(nPDG>0) nPDG=R2;
      else       nPDG=-R2;
#ifdef debug
      G4cout<<"G4QParton::ReDiQADiQ:L1="<<L1<<",R2="<<R2<<",L="<<sPDG<<",R="<<nPDG<<G4endl;
#endif
    }
    else //(L2==R2)
    {
      if(sPDG>0) sPDG=L1;
      else       sPDG=-L1;
      if(nPDG>0) nPDG=R1;
      else       nPDG=-R1;
#ifdef debug
      G4cout<<"G4QParton::ReDiQADiQ:L1="<<L1<<",R1="<<R1<<",L="<<sPDG<<",R="<<nPDG<<G4endl;
#endif
    }
    d1->SetPDGCode(sPDG);             // Reset the left quark
    d2->SetPDGCode(nPDG);            // Reset the right quark
    result=true;
#ifdef debug
    G4cout<<"G4QParton::ReduceDiQADiQ:AfterReduction,L="<<sPDG<<",R="<<nPDG<<G4endl;
#endif
  }
#ifdef debug
  else G4cout<<"-Warning-G4QParton::ReduceDiQADiQ:DQ-aDQ reduction to Q-aQ Failed"<<G4endl;
#endif
  return result;
}
