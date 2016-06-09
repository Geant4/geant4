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
// $Id: G4QString.hh,v 1.2 2006/12/12 11:02:22 mkossov Exp $
// GEANT4 tag $Name: geant4-09-02 $

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

#include "G4ios.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include "G4QHadronVector.hh"
#include "G4QHadronBuilder.hh"
#include "G4QPartonPair.hh"
#include "G4QPartonVector.hh"
#include <algorithm>

class G4QString 
{
 public:

		enum {PROJECTILE  = 1, TARGET  = -1}; // The same as in quark-pair (@@ ? M.K.)

  G4QString(); // formal creation of the string with future excitation
  G4QString(G4QParton* Color,G4QParton* Gluon, G4QParton* AntiColor, G4int Dir=PROJECTILE);
  G4QString(G4QParton* Col,G4QParton* AntiCol, G4int Dir=PROJECTILE);
  G4QString(const G4QString &right);

  ~G4QString();

  // Selectors
  G4int operator==(const G4QString &right) const {return this == &right;}
  G4int operator!=(const G4QString &right) const {return this != &right;}
  const G4ThreeVector& GetPosition() const       {return thePosition;}
  const G4QPartonVector* GetPartonList() const   {return &thePartons;}
  G4QParton* GetGluon() const                    {return thePartons[1];}
  G4QParton* GetGluon(G4int GluonPos) const      {return thePartons[1 + GluonPos];}
  G4QParton* GetLeftParton() const               {return *thePartons.begin();}
  G4QParton* GetRightParton() const              {return *(thePartons.end()-1);}
  G4bool     IsItKinkyString() const             {return (thePartons.size() > 2);}
  G4int      GetDirection() const                {return theDirection;}
  G4bool     IsExcited() const                   {return true;} // Can be absolete? @@ M.K.
  G4LorentzVector Get4Momentum() const;
  //G4QHadron* GetAsQHadron(); //@@ Create a QHadron, corresponding to not excited string
  G4QParton* GetColorParton() const;
  G4QParton* GetAntiColorParton() const;
  G4QParton* GetStableParton() const{ return theStableParton;} // stable at the moment
  G4QParton* GetDecayParton() const{return theDecayParton;}// now involved in fragmentation
  G4int GetDecayDirection() const;
  G4bool FourQuarkString() const; // Checks that the string is qq-(qq)bar string
  G4bool DecayIsQuark() const {return theDecayParton->GetParticleSubType()=="quark";}
  G4bool StableIsQuark() const {return theStableParton->GetParticleSubType()=="quark";}
  G4ThreeVector StablePt(); // Get Pt of the stable quark
  G4ThreeVector DecayPt();  // Get Pt of the decaying quark
  G4double LightConePlus(){return Pplus;}
  G4double LightConeMinus() {return Pminus;}
  G4double LightConeDecay();
  G4LorentzVector GetFragmentation4Mom() const; //!! Instead of Get$Momentum in FragmString
  G4double Mass2() const {return Pplus*Pminus-(Ptleft+Ptright).mag2();}
  G4double Mass() const{return std::sqrt(Mass2());}
  G4bool StopFragmentation() {return FragmentationMass(&G4QHadronBuilder::BuildHighSpin)
                                                 + MassCut > GetFragmentation4Mom().mag();}
  G4bool IsFragmentable() {return FragmentationMass() + MassCut < Mass();}
  G4ThreeVector SampleQuarkPt();
  G4bool SplitLast(G4QString* string, G4QHadronVector* LeftHV, G4QHadronVector* RightHV);

  void Sample4Momentum(G4LorentzVector* Mom, G4double M, G4LorentzVector* aMom,
                       G4double aM, G4double InitialMass);
  // Modifiers
  void SetPosition(const G4ThreeVector& aPosition) {thePosition= aPosition;}
  void SetDirection(G4int dir)                     {theDirection=dir;}
  void LorentzRotate(const G4LorentzRotation& rotation);
  void InsertParton(G4QParton* aParton, const G4QParton* addafter = NULL);
  void Boost(G4ThreeVector& Velocity);
  G4LorentzRotation TransformToCenterOfMass();
  G4LorentzRotation TransformToAlignedCms();
  void DiffString(G4QHadron* aHadron, G4bool isProjectile);
  void ExciteString(G4QParton* Col,G4QParton* AntiCol, G4int Dir);
  G4QHadronVector* FragmentString(G4bool QL); // Fragment String using QGSM=true/LUND=false
  G4QString* CPExcited(); // Creates a string, using only CMS end-partons of the string
  G4QHadronVector* LightFragmentationTest();
  G4double FragmentationMass(G4QHcreate build=0, G4QHadronPair* pdefs=0);
  void CalculateHadronTimePosition(G4double theInitialStringMass, G4QHadronVector*);
  void SetLeftPartonStable();
  void SetRightPartonStable();
  void UpdateString(G4QParton* decay, const G4LorentzVector* mom);
  G4QString(G4QParton* newdecay, const G4LorentzVector* momentum);
  G4QHadron* Splitup(G4bool QL);  // Split Hadron & update the string, QGSM:true,Lund:false
  G4LorentzVector* SplitEandP(G4QHadron* pHadron, G4bool QL); // QGSM:QL=true,Lund:QL=false
  G4QParton* CreateParton(G4int PDGcode)                   {return new G4QParton(PDGcode);}
  G4QPartonPair CreatePartonPair(G4int NeedParticle, G4bool AllowDiquarks=true);
  G4QHadron* QuarkSplitup(G4QParton* decay, G4QParton* &created);
  G4QHadron* DiQuarkSplitup(G4QParton* decay, G4QParton* &created);
  G4int SampleQuarkFlavor() {return (1+G4int(G4UniformRand()/StrangeSuppress));} // ? M.K.

  // Static functions
  static void SetParameters(G4double mCut, G4double clustM, G4double sigQT, G4double DQSup,
   G4double DQBU, G4double smPar, G4double SSup, G4double SigPt, G4int SLmax, G4int CLmax);

 private:
  // Private functions
  G4ThreeVector GaussianPt(G4double widthSquare, G4double maxPtSquare) const;
  G4double GetLundLightConeZ(G4double zmin, G4double zmax, G4int PartonEncoding,
                             G4QHadron* pHadron, G4double Px, G4double Py);
  G4double GetQGSMLightConeZ(G4double zmin, G4double zmax, G4int PartonEncoding,
                             G4QHadron* pHadron, G4double Px, G4double Py);

  // Static parameters
  // Parameters of Longitudinal String Fragmentation
  static G4double MassCut;           // minimum mass cut for the string
  static G4double ClusterMass;       // minimum cluster mass for the fragmentation
  static G4double SigmaQT;           // quark transverse momentum distribution parameter 
  static G4double DiquarkSuppress;   // is Diquark suppression parameter  
  static G4double DiquarkBreakProb;  // is Diquark breaking probability 
  static G4double SmoothParam;       // QGS model parameter
  static G4double StrangeSuppress;   // Strangeness suppression parameter
	 static G4double widthOfPtSquare;	  // width^2 of pt for string excitation
  static G4int StringLoopInterrupt;  // String fragmentation LOOP limit 
  static G4int ClusterLoopInterrupt; // Cluster fragmentation LOOP limit 

  // Body
  G4int         theDirection;  // must be 1 or -1 (PROJECTILE or TARGET)
  G4ThreeVector thePosition;   // Defined by the first quark position
  G4QPartonVector thePartons;  // would like initial capacity for 3 Partons (? M.K.)
  G4QHadronBuilder* hadronizer;// Hadronizer of hodrons out of partons
  G4ThreeVector Ptleft,Ptright;// Pt (px,py) for partons (pz ignored!)
  G4double Pplus, Pminus;      // p-, p+ of string, Plus is assigned to Left!
  G4QParton* theStableParton;  // Parton on the stable side of the string
  G4QParton* theDecayParton;   // Parton on the decaying part of the string
  enum DecaySide {None, Left, Right}; // @@ To have two         @@ Leav   :  1=Left
  DecaySide decaying;                 // @@   it's too much     @@  only  :  0=Unknown
  G4int     SideOfDecay;              // @@     of a good thing @@   one! : -1=Right
};

#endif
