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
// $Id: G4Parton.hh,v 1.14 2002/12/12 19:17:57 gunter Exp $
// GEANT4 tag $Name: geant4-05-01 $
//

#ifndef G4Parton_h
#define G4Parton_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4Parton ----------------
//             by Gunter Folger, June 1998.
//       class for Parton (inside a string) used by Parton String Models
// ------------------------------------------------------------

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "g4std/iostream"
#include "G4ParticleTable.hh"
#include "Randomize.hh"

class G4Parton
{
   public:
      G4Parton()
      {
        // CAUTION: 
        // this is a preliminary definition yielding u and d quarks only!
        //
        PDGencoding=(G4int)(2.*G4UniformRand()); 
        theColour = (G4int)(3.*G4UniformRand())+1;
        theIsoSpinZ = ((G4int)(G4UniformRand()))-0.5;
        theSpinZ = ((G4int)(G4UniformRand()))-0.5;
      }
      
      G4Parton(G4int PDGencoding);
      G4Parton(const G4Parton &right);

      ~G4Parton();

      const G4Parton & operator=(const G4Parton &right);

      int operator==(const G4Parton &right) const;

      int operator!=(const G4Parton &right) const;
      
      G4int GetPDGcode() const;
      
      G4ParticleDefinition * GetDefinition();

      void DefineMomentumInZ(G4double aLightConeMomentum, G4bool aDirection);      
      
      const G4ThreeVector & GetPosition()const;
      void SetPosition(const G4ThreeVector &aPosition);

      const G4LorentzVector & Get4Momentum() const; 
      void Set4Momentum(const G4LorentzVector & aMomentum);
      
      void SetX(G4double anX) { theX = anX; }
      
      void SetColour(G4int aColour) {theColour = aColour;}
      G4int GetColour() {return theColour;}
      
      void SetIsoSpinZ(G4double anIsoSpinZ) {theIsoSpinZ = anIsoSpinZ;}
      G4double GetIsoSpinZ() {return theIsoSpinZ;}

      void SetSpinZ(G4double aSpinZ) {theSpinZ = aSpinZ;}
      G4double GetSpinZ() {return theSpinZ;}
      
   private:
      G4double GetMass();
      
    private:  
      G4int PDGencoding;
      G4ParticleDefinition * theDefinition;
      G4LorentzVector theMomentum;
      G4ThreeVector   thePosition;
      
      G4int theColour;
      G4double theIsoSpinZ;
      G4double theSpinZ;
      
      G4double theX;
      
};

inline int G4Parton::operator==(const G4Parton &right) const
{
	return this==&right;
}	

inline int G4Parton::operator!=(const G4Parton &right) const
{
	return this!=&right;
}

inline G4int G4Parton::GetPDGcode() const
{
	return PDGencoding;
}
	
inline const G4ThreeVector & G4Parton::GetPosition() const
{
	return thePosition;
}

inline void G4Parton::SetPosition(const G4ThreeVector &aPosition)
{
	thePosition=aPosition;
}


inline const G4LorentzVector & G4Parton::Get4Momentum() const
{
	return theMomentum;
}

inline void G4Parton::Set4Momentum(const G4LorentzVector & aMomentum)
{
	theMomentum=aMomentum;
}


inline
G4double G4Parton::GetMass()
{
	return theDefinition->GetPDGMass();
}

inline
G4ParticleDefinition * G4Parton::GetDefinition()
{
	return theDefinition;
}


#endif
