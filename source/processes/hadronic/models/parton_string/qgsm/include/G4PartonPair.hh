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
    enum {
	   PROJECTILE = 1,
	   TARGET = -1
	  };
  public:
    G4PartonPair(G4Parton* P1, G4Parton* P2, G4int Type, G4int Direction);
    ~G4PartonPair();

  private:
    G4PartonPair(const G4PartonPair &right);
    int operator==(const G4PartonPair &right) const;
    int operator!=(const G4PartonPair &right) const;

  public:
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

