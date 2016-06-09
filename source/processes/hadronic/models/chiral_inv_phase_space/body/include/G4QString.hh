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
// $Id$

#ifndef G4QString_h
#define G4QString_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4QString ----------------
//       by Mikhail Kosov, October 2006.
//       class for an excited string used by Parton String Models
//   For comparison mirror member functions are taken from G4 classes:
//   G4FragmentingString
//   G4ExcitedStringDecay
//
// ------------------------------------------------------------
// Short description: If partons from the G4QPartonPair are close in
// rapidity, they create Quasmons, but if they are far in the rapidity
// space, they can not interact directly. Say the bottom parton (quark)
// has rapidity 0, and the top parton (antiquark) has rapidity 8, then
// the top quark splits in two by radiating gluon, and each part has
// rapidity 4, then the gluon splits in quark-antiquark pair (rapidity
// 2 each), and then the quark gadiates anothe gluon and reachs rapidity
// 1. Now it can interact with the bottom antiquark, creating a Quasmon
// or a hadron. The intermediate partons is the string ladder.
// ---------------------------------------------------------------------

#include "G4ios.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include "G4QHadronVector.hh"
#include "G4QPartonPair.hh"
#include "G4QPartonVector.hh"
#include <algorithm>

class G4QString 
{
 public:

  enum {PROJECTILE  = 1, TARGET  = -1}; // The same as in quark-pair (@@ ? M.K.)

  G4QString(); // formal creation of the string with future excitation
  G4QString(G4QParton* Color,G4QParton* Gluon, G4QParton* AntiColor, G4int Dir=PROJECTILE);
  G4QString(G4QParton* Col, G4QParton* AntiCol, G4int Dir=PROJECTILE);
  G4QString(G4QPartonPair* ColAntiCol);
  G4QString(const G4QString &right);
  G4QString(const G4QString &old, G4QParton* newdecay, const G4LorentzVector *momentum);
  G4QString(const G4QString &old, G4QParton* newdecay);
  G4QString(G4QParton* newdecay, const G4LorentzVector* momentum);

  ~G4QString();

  // Selectors
  G4int operator==(const G4QString &right) const {return this == &right;}
  G4int operator!=(const G4QString &right) const {return this != &right;}
  const G4ThreeVector& GetPosition() const       {return thePosition;}
  const G4QPartonVector* GetPartonList() const   {return &thePartons;}
  G4QParton* GetLeftParton() const               {return *thePartons.begin();}
  G4QParton* GetRightParton() const              {return *(thePartons.end()-1);}
  G4int      GetDirection() const                {return theDirection;}
  G4LorentzVector Get4Momentum() const;
  G4QContent GetQC() const {return GetLeftParton()->GetQC()+GetRightParton()->GetQC();}
  G4int      GetCharge() const {return GetQC().GetCharge();}
  G4int      GetBaryonNumber() const {return GetQC().GetBaryonNumber();}
  G4int      GetStrangeness() const {return GetQC().GetStrangeness();}
  G4int GetDecayDirection() const;
  G4bool DecayIsQuark() const {return theDecayParton->GetType()==1;}
  G4bool StableIsQuark() const {return theStableParton->GetType()==1;}
  G4ThreeVector DecayPt();  // Get Pt of the decaying quark @@ Called once
  G4double Mass2() const { return Pplus*Pminus-(Ptleft+Ptright).mag2();}
  G4double Mass() const  // @@ Very dangerous! USE ONLY FORE THE LIGHT CONE ALGORITHM !!
  {
    G4double  mass2=Mass2();
    if(mass2>0) return std::sqrt(Mass2());
    else     return 0.; // @@ Make Warning
  }
  G4bool StopFragmentation()
  {
    G4LorentzVector mom(Ptleft+Ptright, 0.5*(Pplus+Pminus));
    mom.setPz(0.5*(Pplus-Pminus));
    return FragmentationMass(1) + MassCut > mom.mag();
  }
  G4bool IsFragmentable() {return FragmentationMass() + MassCut < Mass();} // @@ Mass() ?
  G4ThreeVector SampleQuarkPt(); // @@ Called once
  G4QHadron* CreateHadron(G4QParton* black, G4QParton* white);
  G4QHadron* CreateLowSpinHadron(G4QParton* black, G4QParton* white);
  G4QHadron* CreateHighSpinHadron(G4QParton* black, G4QParton* white);

  // Modifiers
  void SetPosition(const G4ThreeVector& aPosition){thePosition= aPosition;}
  void SetDirection(G4int dir)       {if(dir==1 || dir==-1) theDirection=dir;}
  void SetLeftParton(G4QParton* LP)  {thePartons[0]=LP;} // !! Not deleting the substituty
  void SetRightParton(G4QParton* RP) {thePartons.pop_back(); thePartons.push_back(RP);}
  void KillString()                  {theDirection=0;} // @@ Can be absolete
  void LorentzRotate(const G4LorentzRotation& rotation);
  //void InsertParton(G4QParton* aParton, const G4QParton* addafter = NULL);
  void Boost(G4ThreeVector& Velocity);
  G4LorentzRotation TransformToAlignedCms(); // @@ caled once
  //void DiffString(G4QHadron* aHadron, G4bool isProjectile); @@ Mast be used!!
  void ExciteString(G4QParton* Col,G4QParton* AntiCol, G4int Dir);
  G4QHadronVector* FragmentString(G4bool QL); // Fragment String using QGSM=true/LUND=false
  G4QHadronVector* LightFragmentationTest();
  G4double FragmentationMass(G4int HighSpin = 0, G4QHadronPair* pdefs = 0);
  void SetLeftPartonStable();
  void SetRightPartonStable();
  G4QHadron* Splitup(G4bool QL);
  G4LorentzVector* SplitEandP(G4QHadron* pHadron, G4bool QL); // QGSM:QL=true,Lund:QL=false
  G4QPartonPair CreatePartonPair(G4int NeedParticle, G4bool AllowDiquarks=true);
  G4QHadron* QuarkSplitup(G4QParton* decay, G4QParton* &created);
  G4QHadron* DiQuarkSplitup(G4QParton* decay, G4QParton* &created);
  G4int SampleQuarkFlavor() {return (1+G4int(G4UniformRand()/StrangeSuppress));} // ? M.K.

  // Static functions
  static void SetParameters(G4double mCut, G4double sigQT, G4double DQSup, G4double DQBU,
                            G4double smPar, G4double SSup, G4double SigPt);

 private:
  enum Spin {SpinZero=1, SpinHalf=2, SpinOne=3, SpinThreeHalf=4};
  // Private functions
  G4QHadron* CreateMeson(G4QParton* black, G4QParton* white, Spin spin);
  G4QHadron* CreateBaryon(G4QParton* black,G4QParton* white, Spin spin);
  G4ThreeVector GaussianPt(G4double widthSquare, G4double maxPtSquare) const;
  G4double GetLundLightConeZ(G4double zmin, G4double zmax, G4int PartonEncoding,
                             G4QHadron* pHadron, G4double Px, G4double Py);
  G4double GetQGSMLightConeZ(G4double zmin, G4double zmax, G4int PartonEncoding,
                             G4QHadron* pHadron, G4double Px, G4double Py);

  // Static parameters
  // Parameters of Longitudinal String Fragmentation
  static G4double MassCut;           // minimum mass cut for the string
  static G4double SigmaQT;           // quark transverse momentum distribution parameter 
  static G4double DiquarkSuppress;   // is Diquark suppression parameter  
  static G4double DiquarkBreakProb;  // is Diquark breaking probability 
  static G4double SmoothParam;       // QGS model parameter
  static G4double StrangeSuppress;   // Strangeness suppression parameter
  static G4double widthOfPtSquare;   // width^2 of pt for string excitation

  // Body
  G4int         theDirection;        // must be 1 (PROJECTILE) or -1 (TARGET), 0 - DEAD
  G4ThreeVector thePosition;         // Defined by the first quark position
  G4QPartonVector thePartons;        // Partons on the ends of the string @@ Use PartonPair
  G4ThreeVector Ptleft,Ptright;      // Pt (px,py) for partons (pz ignored!)
  G4double Pplus, Pminus;            // p-, p+ of string, Plus is assigned to Left!
  G4QParton* theStableParton;        // Parton on the stable side of the string
  G4QParton* theDecayParton;         // Parton on the decaying part of the string
  enum DecaySide {None, Left, Right};// @@ To have two         @@ Leav   :  1=Left
  DecaySide decaying;                // @@   it's too much     @@  only  :  0=Unknown
  G4int     SideOfDecay;             // @@     of a good thing @@   one! : -1=Right
};

#endif
