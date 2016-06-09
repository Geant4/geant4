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
#ifndef G4PartonPair_h
#define G4PartonPair_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4Parton.hh"
#include "G4PartonVector.hh"

class G4PartonPair 
{
  public:
      enum {
           DIFFRACTIVE = 1,
           SOFT = 2,
           HARD  = 3 
           };
      enum
           {
           PROJECTILE = 1,
           TARGET = -1
           };   
  public:
      G4PartonPair(G4Parton* P1, G4Parton* P2, G4int Type, G4int Direction);
      G4PartonPair(const G4PartonPair &right);
      ~G4PartonPair();
      
      int operator==(const G4PartonPair &right) const;
      int operator!=(const G4PartonPair &right) const;

      void  SetPartons(G4Parton* P1, G4Parton* P2);
      void  SetCollisionType(G4int Type);
      G4int GetCollisionType();
      G4Parton* GetParton1(void);
      G4Parton* GetParton2(void);
      G4int GetDirection();
      
      
  private:
      G4Parton* Parton1;  
      G4Parton* Parton2;  
      G4int     CollisionType;
      G4int     Direction;

};

inline int G4PartonPair::operator==(const G4PartonPair &right) const
    {
    return (CollisionType == right.CollisionType && 
       *Parton1 == *right.Parton1 &&
       *Parton2 == *right.Parton2)? 1: 0;
    }

inline int G4PartonPair::operator!=(const G4PartonPair &right) const
    {
    return (CollisionType == right.CollisionType && 
       *Parton1 == *right.Parton1 &&
       *Parton2 == *right.Parton2)? 0: 1;
    }

inline G4Parton* G4PartonPair::GetParton1(void)
    {
    return Parton1;
    }

inline G4Parton* G4PartonPair::GetParton2(void)
    {
    return Parton2;
    }

inline void G4PartonPair::SetCollisionType(G4int Type)
    {
    CollisionType = Type;
    }

inline G4int G4PartonPair::GetCollisionType()
    {
    return CollisionType;
    }
    
inline G4int G4PartonPair::GetDirection()
    {
    return Direction;
    }


#endif


