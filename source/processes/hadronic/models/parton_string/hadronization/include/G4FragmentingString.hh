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
// $Id: G4FragmentingString.hh 102717 2017-02-20 10:37:13Z gcosmo $
//

#ifndef G4FragmentingString_h
#define G4FragmentingString_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4FragmentingString ----------------
//             by Gunter Folger, September 2001.
//       class for an excited string used in Fragmention
// ------------------------------------------------------------

#include "G4ios.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include "G4ParticleDefinition.hh"

class G4ExcitedString;

class G4FragmentingString 
{

  public:
     
      G4FragmentingString(const G4FragmentingString &right);
      G4FragmentingString(const G4ExcitedString &excited);
      G4FragmentingString(const G4FragmentingString &old,
			  G4ParticleDefinition * newdecay,
			  const G4LorentzVector *momentum);
      G4FragmentingString(const G4FragmentingString &old,
			  G4ParticleDefinition * newdecay);
			  
      ~G4FragmentingString();

      G4FragmentingString& operator=(const G4FragmentingString &);
      int operator==(const G4FragmentingString &right) const;

      int operator!=(const G4FragmentingString &right) const;

      
      G4LorentzVector Get4Momentum() const;

      G4ThreeVector StablePt();
      G4ThreeVector DecayPt();
      
      G4double LightConePlus();
      G4double LightConeMinus();
      G4double LightConeDecay();
   
      G4double Mass() const;
      G4double Mass2() const;
      G4double MassT2() const;
      
      G4ParticleDefinition* GetLeftParton(void) const;
      G4ParticleDefinition* GetRightParton(void) const;
      
      G4ParticleDefinition* GetStableParton() const; // stable at the moment
      G4ParticleDefinition* GetDecayParton() const;  // currently involved in fragmentation

      void SetLeftPartonStable();
      void SetRightPartonStable();

      G4int GetDecayDirection() const;
            
      G4bool    DecayIsQuark();
      G4bool    StableIsQuark();
      G4bool    FourQuarkString(void) const;

G4LorentzVector   GetPstring();
G4LorentzVector   GetPleft();
void              SetPleft(G4LorentzVector a4momentum);
G4LorentzVector   GetPright();
void              SetPright(G4LorentzVector a4momentum);
void              LorentzRotate(const G4LorentzRotation & rotation);
G4LorentzRotation TransformToCenterOfMass();
G4LorentzRotation TransformToAlignedCms();
void              Boost(G4ThreeVector& Velocity);


  private:

      G4ParticleDefinition *LeftParton, *RightParton; 
      G4ThreeVector Ptleft,Ptright;    // Pt (px,py) for partons (pz ignored!)
      G4double Pplus, Pminus;        // p-, p+ of string, Plus ass. to Left!
  
      G4ParticleDefinition * theStableParton, * theDecayParton;
      
      G4LorentzVector Pstring, Pleft, Pright;
      enum DecaySide { None, Left, Right };
      DecaySide decaying; 
};

inline
int G4FragmentingString::operator==(const G4FragmentingString &right) const
{
	return this == &right;
}

inline
int G4FragmentingString::operator!=(const G4FragmentingString &right) const
{
	return this != &right;
}


inline
G4ParticleDefinition * G4FragmentingString::GetStableParton() const
{
     return  theStableParton;
}	

inline
G4ParticleDefinition * G4FragmentingString::GetDecayParton() const
{
     return  theDecayParton;
}	

inline
G4ParticleDefinition* G4FragmentingString::GetLeftParton(void) const
    {
    return LeftParton; 
    }

inline
G4ParticleDefinition* G4FragmentingString::GetRightParton(void) const
    {
    return RightParton; 
    }

//+++++++++++++++++++++++++++
inline
void G4FragmentingString::LorentzRotate(const G4LorentzRotation & rotation)
{
     SetPleft(rotation*Pleft);
     SetPright(rotation*Pright);
     Pstring = Pleft+Pright;
Ptleft =Pleft.vect();  Ptleft.setZ(0.);
Ptright=Pright.vect(); Ptright.setZ(0.);
Pplus =Pstring.plus();
Pminus=Pstring.minus();
}

inline
G4LorentzRotation G4FragmentingString::TransformToCenterOfMass()
{
     G4LorentzVector momentum=Pstring;
     G4LorentzRotation toCMS(-1*momentum.boostVector());

     Pleft   *= toCMS;
     Pright  *= toCMS;
     Pstring *= toCMS;
Ptleft =Pleft.vect();  Ptleft.setZ(0.);
Ptright=Pright.vect(); Ptright.setZ(0.);
Pplus =Pstring.plus();
Pminus=Pstring.minus();
     return toCMS;
}

inline
G4LorentzRotation G4FragmentingString::TransformToAlignedCms()
{
     G4LorentzVector momentum=Pstring;
     G4LorentzRotation toAlignedCms(-1*momentum.boostVector());

     momentum= toAlignedCms* Pleft;
     toAlignedCms.rotateZ(-1*momentum.phi());
     toAlignedCms.rotateY(-1*momentum.theta());

     Pleft   *= toAlignedCms;
     Pright  *= toAlignedCms;
     Pstring *= toAlignedCms;

Ptleft  = G4ThreeVector(0.,0.,0.);
Ptright = G4ThreeVector(0.,0.,0.);
Pplus  = Pstring.plus();
Pminus = Pstring.minus();

     return toAlignedCms;
}

inline
void G4FragmentingString::SetPleft(G4LorentzVector a4momentum)
{    Pleft = a4momentum;}

inline
void G4FragmentingString::SetPright(G4LorentzVector a4momentum)
{    Pright = a4momentum;}
#endif


