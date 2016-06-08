// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ExcitedString.hh,v 1.5 1999/12/15 14:52:50 gunter Exp $
// GEANT4 tag $Name: geant4-03-00 $
//

#ifndef G4ExcitedString_h
#define G4ExcitedString_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
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
     
      G4ExcitedString(const G4ExcitedString &right);

      ~G4ExcitedString();

      int operator==(const G4ExcitedString &right) const;

      int operator!=(const G4ExcitedString &right) const;

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
	for ( G4int index=0; index < thePartons.entries() ; index++ )
	{
	    momentum += thePartons[index]->Get4Momentum();
	}
	return momentum;
}

inline
void G4ExcitedString::LorentzRotate(const G4LorentzRotation & rotation)
{
	for ( G4int index=0; index < thePartons.entries() ; index++ )
	{
	    thePartons[index]->Set4Momentum(rotation*thePartons[index]->Get4Momentum());
	}
}

inline
void G4ExcitedString::InsertParton(G4Parton *aParton, const G4Parton * addafter)
{

	G4int insert_index=1;
	
	if ( addafter != NULL ) 
	{
	   insert_index=thePartons.index(addafter);
	   if (insert_index == G4std::string::npos)		// No object addafter in thePartons
	   {
	   	G4Exception("G4ExcitedString::InsertParton called with invalid second argument");
	   }
	}
	
	thePartons.insertAt(insert_index+1, aParton);
} 

inline
G4LorentzRotation G4ExcitedString::TransformToCenterOfMass()
{
	G4LorentzVector momentum=Get4Momentum();
	G4LorentzRotation toCms(-1*momentum.boostVector());

	for ( G4int index=0; index < thePartons.entries() ; index++ )
	{
	    momentum=toCms * thePartons[index]->Get4Momentum();
	    thePartons[index]->Set4Momentum(momentum);
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
	
	for ( G4int index=0; index < thePartons.entries() ; index++ )
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


