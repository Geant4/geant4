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
// $Id: G4FragmentingString.hh 106967 2017-10-31 08:41:49Z gcosmo $
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

