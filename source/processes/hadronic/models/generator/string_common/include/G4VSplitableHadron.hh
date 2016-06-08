// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VSplitableHadron.hh,v 1.3 1999/12/15 14:52:46 gunter Exp $
// GEANT4 tag $Name: geant4-03-01 $
//

#ifndef G4VSplitableHadron_h
#define G4VSplitableHadron_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4VSplitableHadron----------------
//             by Gunter Folger, June 1998.
//       class storing an interacting particle. Used by Parton String Models.
// ------------------------------------------------------------

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ReactionProduct.hh"
class G4Nucleon;
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
class G4Parton;
class G4VKineticNucleon;

#include "g4rw/tpordvec.h"

class G4VSplitableHadron 
{

  public:
      G4VSplitableHadron();
      G4VSplitableHadron(const G4ReactionProduct & aPrimary);
      G4VSplitableHadron(const G4Nucleon & aNucleon);
      G4VSplitableHadron(const G4VKineticNucleon * aNucleon);

      virtual ~G4VSplitableHadron();


      int operator==(const G4VSplitableHadron &right) const;
      int operator!=(const G4VSplitableHadron &right) const;

      void Set4Momentum(const G4LorentzVector &a4Momentum);
      const G4LorentzVector & Get4Momentum() const;

      void SetDefinition(G4ParticleDefinition *aDefinition);
      G4ParticleDefinition * GetDefinition() const;
      
      void IncrementCollisionCount(G4int aCount);

      void SetPosition(const G4ThreeVector &aPosition);
      const G4ThreeVector & GetPosition() const;

      virtual void SplitUp() = 0;
      virtual G4Parton * GetNextParton() = 0 ;
      virtual G4Parton * GetNextAntiParton() = 0 ;


  protected:
      G4int GetSoftCollisionCount();
      
  private:

      G4VSplitableHadron(const G4VSplitableHadron &right);
      const G4VSplitableHadron & operator=(const G4VSplitableHadron &right);

  private:

      G4ParticleDefinition *theDefinition;

      G4LorentzVector the4Momentum;

      G4ThreeVector thePosition;
      G4int theCollisionCount;

};

inline G4int G4VSplitableHadron::GetSoftCollisionCount()
{
return theCollisionCount;
}

inline void G4VSplitableHadron::Set4Momentum(const G4LorentzVector &a4Momentum)
{
	the4Momentum=a4Momentum;
}

inline const G4LorentzVector & G4VSplitableHadron::Get4Momentum() const
{
	return the4Momentum;
}

inline void G4VSplitableHadron::SetDefinition(G4ParticleDefinition *aDefinition)
{
	theDefinition=aDefinition;
}

inline G4ParticleDefinition * G4VSplitableHadron::GetDefinition() const
{
	return theDefinition;
}

inline void G4VSplitableHadron::IncrementCollisionCount(G4int aCount)
{
	theCollisionCount += aCount;
}


inline void G4VSplitableHadron::SetPosition(const G4ThreeVector &aPosition)
{
	thePosition=aPosition;
}

inline const G4ThreeVector & G4VSplitableHadron::GetPosition() const
{
	return thePosition;
}



#endif


