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
#ifndef G4QPartonPair_h
#define G4QPartonPair_h 1
//
// $Id: G4QPartonPair.hh,v 1.2 2009-02-23 09:49:24 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4QParton ----------------
//             by Mikhail Kossov, Oct 2006.
// class for PartonPair (hadron) used by Parton String Models
// ------------------------------------------------------------
// Short description: Each Quasmon String has a pair of partons
// (quark/diquark-partons) on its ends. During the hadronization
// procedure the rapidity gap between partons shrinks, but the
// parton pair still exists, while it is converted to the final
// meson (quaek-antiquark) or baryon (quark-diquark).
// --------------------------------------------------------------
//
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4QParton.hh"

class G4QPartonPair 
{
 public:
  enum { DIFFRACTIVE = 1, SOFT = 2, HARD  = 3};
  enum { PROJECTILE = 1, TARGET = -1};
 public:
  G4QPartonPair(G4QParton* P1, G4QParton* P2, G4int Type=0, G4int Direction=0);
  ~G4QPartonPair();
  G4int operator==(const G4QPartonPair &right) const
  {
    return (CollisionType == right.CollisionType &&
            *Parton1 == *right.Parton1 && *Parton2 == *right.Parton2) ? 1: 0;
  }
  G4int operator!=(const G4QPartonPair &right) const
  {
    return (CollisionType == right.CollisionType &&
            *Parton1 == *right.Parton1 && *Parton2 == *right.Parton2) ? 0: 1;
  }
  // Modifiers
  void  SetPartons(G4QParton* P1, G4QParton* P2) {Parton1=P1; Parton2=P2;}
  void  SetCollisionType(G4int Type)             {CollisionType = Type;}
  // Selectors
  G4int      GetCollisionType()                  {return CollisionType;}
  G4QParton* GetParton1()                        {return Parton1;}
  G4QParton* GetParton2()                        {return Parton2;}
  G4int      GetDirection()                      {return Direction;}
      
 private:
  // Body
  G4QParton* Parton1;  
  G4QParton* Parton2;  
  G4int      CollisionType;
  G4int      Direction;
};

#endif
