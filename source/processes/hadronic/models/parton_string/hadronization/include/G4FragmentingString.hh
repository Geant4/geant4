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
// $Id: G4FragmentingString.hh,v 1.1 2003/10/07 11:25:40 hpw Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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

      ~G4FragmentingString();

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


  private:

      G4ParticleDefinition *LeftParton, *RightParton; 
      G4ThreeVector Ptleft,Ptright;    // Pt (px,py) for partons (pz ignored!)
      G4double Pplus, Pminus;        // p-, p+ of string, Plus ass. to Left!
  
      G4ParticleDefinition * theStableParton, * theDecayParton;
      
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



#endif


