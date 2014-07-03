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

#ifndef G4ExcitedString_h
#define G4ExcitedString_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4ExcitedString ----------------
//             by Gunter Folger, June 1998.
//       class for an excited string used by Parton String Models
// ------------------------------------------------------------

#include "G4ios.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include "G4Parton.hh"
#include "G4PartonVector.hh"
#include "G4KineticTrack.hh"
#include "G4HadronicException.hh"
#include <algorithm>

class G4ExcitedString 
{

  public:

      enum {
           PROJECTILE  = 1, 
           TARGET  = -1 
           };

      G4ExcitedString(G4Parton* Color, G4Parton* Gluon,  G4Parton* AntiColor, G4int Direction=PROJECTILE);
      G4ExcitedString(G4Parton* Color, G4Parton* AntiColor, G4int Direction=PROJECTILE);
      G4ExcitedString(G4KineticTrack * atrack);
     

      ~G4ExcitedString();

private:
      G4ExcitedString(const G4ExcitedString &right);
      G4ExcitedString& operator= (const G4ExcitedString &right);

public:
      int operator==(const G4ExcitedString &right) const;
      int operator!=(const G4ExcitedString &right) const;

public:
      G4double GetTimeOfCreation() const;

      void  SetTimeOfCreation(G4double aTime);

      const G4ThreeVector & GetPosition() const;

      void SetPosition(const G4ThreeVector &aPosition); 

      const G4PartonVector * GetPartonList() const;

      G4LorentzVector Get4Momentum() const;
      void LorentzRotate(const G4LorentzRotation & rotation);

      void InsertParton(G4Parton * aParton, const G4Parton * addafter = NULL);

      G4LorentzRotation TransformToCenterOfMass();
      G4LorentzRotation TransformToAlignedCms();


      void      Boost(G4ThreeVector& Velocity);

      G4Parton* GetColorParton(void) const;
      G4Parton* GetGluon(void) const;
      G4Parton* GetAntiColorParton(void) const;
      G4Parton* GetGluon(G4int GluonPos) const;
      
      G4KineticTrack * GetKineticTrack() const;

      G4Parton* GetLeftParton(void) const;
      G4Parton* GetRightParton(void) const;

      G4bool    IsItKinkyString(void) const;
      G4int     GetDirection(void) const;
      
      G4bool    IsExcited() const;


  private:

      G4int    theDirection;  // must be 1 or -1 (PROJECTILE or TARGET)
      G4double theTimeOfCreation;
      G4ThreeVector thePosition;
      G4PartonVector thePartons;  // would like initial capacity for 3 Partons.
      G4KineticTrack* theTrack;

};

inline
int G4ExcitedString::operator==(const G4ExcitedString &right) const
{
	return this == &right;
}

inline
int G4ExcitedString::operator!=(const G4ExcitedString &right) const
{
	return this != &right;
}

inline
G4double G4ExcitedString::GetTimeOfCreation() const
{
	return theTimeOfCreation;
}

inline
void G4ExcitedString::SetTimeOfCreation(G4double aTime)
{
	theTimeOfCreation=aTime;
}

inline
const G4ThreeVector & G4ExcitedString::GetPosition() const 
{
	return thePosition;
}

inline
void G4ExcitedString::SetPosition(const G4ThreeVector &aPosition)
{
	thePosition= aPosition;
}

inline
G4LorentzVector G4ExcitedString::Get4Momentum() const
{
	G4LorentzVector momentum;
        if ( IsExcited() )
        {
	    for ( unsigned int index=0; index < thePartons.size() ; index++ )
	    {
	        // std::cout << "HPW "<<thePartons[index]->Get4Momentum()<<std::endl;
	        momentum += thePartons[index]->Get4Momentum();
	    }
        }
        else
        {
            momentum=theTrack->Get4Momentum();
        }
	return momentum;
}

inline
void G4ExcitedString::LorentzRotate(const G4LorentzRotation & rotation)
{
        if ( IsExcited() )
        {
	    for ( unsigned int index=0; index < thePartons.size() ; index++ )
	    {
	        thePartons[index]->Set4Momentum(rotation*thePartons[index]->Get4Momentum());
	    }
        }
        else
        {
            theTrack->Set4Momentum(rotation*theTrack->Get4Momentum());
        }
}

inline
void G4ExcitedString::InsertParton(G4Parton *aParton, const G4Parton * addafter)
{

	G4PartonVector::iterator insert_index;
	
	if ( addafter != NULL ) 
	{
	    insert_index=std::find(thePartons.begin(), thePartons.end(), addafter);
	    if ( insert_index == thePartons.end() )		// No object addafter in thePartons
	    {
	        G4String text = "G4ExcitedString::InsertParton called with invalid second argument";
  	        throw G4HadronicException(__FILE__, __LINE__, text);
	    }
	}
	
	thePartons.insert(insert_index+1, aParton);
} 

inline
G4LorentzRotation G4ExcitedString::TransformToCenterOfMass()
{
	G4LorentzVector momentum=Get4Momentum();
	G4LorentzRotation toCms(-1*momentum.boostVector());

        if ( IsExcited() )
        {
	    for ( unsigned int index=0; index < thePartons.size() ; index++ )
	    {
	        momentum=toCms * thePartons[index]->Get4Momentum();
	        thePartons[index]->Set4Momentum(momentum);
	    }
        }
        else
        {
            momentum*=toCms;
            theTrack->Set4Momentum(momentum);
        }
	return toCms;
}

inline
G4LorentzRotation G4ExcitedString::TransformToAlignedCms()
{
	G4LorentzVector momentum=Get4Momentum();
	G4LorentzRotation toAlignedCms(-1*momentum.boostVector());

	momentum= toAlignedCms* thePartons[0]->Get4Momentum();
	toAlignedCms.rotateZ(-1*momentum.phi());
	toAlignedCms.rotateY(-1*momentum.theta());
	
	for ( unsigned int index=0; index < thePartons.size() ; index++ )
	{
	    momentum=toAlignedCms * thePartons[index]->Get4Momentum();
	    thePartons[index]->Set4Momentum(momentum);
	}
	return toAlignedCms;
}


inline 
const G4PartonVector * G4ExcitedString::GetPartonList() const
{
	return &thePartons;
}

inline 
G4KineticTrack * G4ExcitedString::GetKineticTrack() const
{
	return theTrack;
}

inline 
G4bool G4ExcitedString::IsExcited() const
{
	return theTrack == 0;
}


#endif
