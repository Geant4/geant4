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
// $Id: G4VHit.hh,v 1.6 2001/07/13 15:00:17 gcosmo Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//

#ifndef G4VHit_h
#define G4VHit_h 1

#include "globals.hh"

// class description:
//
//  This is the base class of hit object. The user should derive this
// base class to make his/her own hit class. Two virtual method Draw()
// and Print() can be implemented if the user wants these functionarities.
//  If a concrete hit class is used as a transient class, G4Allocator
// must be used.

class G4VHit 
{

  public:
      G4VHit();
      virtual ~G4VHit();

      G4int operator==(const G4VHit &right) const;

      virtual void Draw();
      virtual void Print();

};

#endif

