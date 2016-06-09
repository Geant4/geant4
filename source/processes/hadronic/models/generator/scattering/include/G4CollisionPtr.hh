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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name:     G4CollisionPtr
//
//      Author:        Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#ifndef G4COLLISIONPTR_HH
#define G4COLLISIONPTR_HH

#include "globals.hh"

class G4VCollision;

class G4CollisionPtr
{
public:

  // This is a wrapper class to store pointers to G4VCollision in vectors
  // Constructor 
  G4CollisionPtr(G4VCollision* x = 0);

  //Destructor
  ~G4CollisionPtr() { }

  // Copy constructor
  G4CollisionPtr(const G4CollisionPtr& xw) : x_(xw.x_) { }

  // Operators

  const G4VCollision* operator() () const;
  G4VCollision* operator() ();

  G4CollisionPtr& operator= (const G4CollisionPtr& xw);

  G4bool operator== (const G4CollisionPtr& right) const;

  G4bool operator< (const G4CollisionPtr& ) { return false; }  

private:  

  G4VCollision* x_;

};

#endif
