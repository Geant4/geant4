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
// $Id: G4VDigi.hh,v 1.7 2001/07/13 15:00:12 gcosmo Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//

#ifndef G4VDigi_h
#define G4VDigi_h 1

#include "globals.hh"

// class description:
//
//  This is the base class of digi object. The user should derive this
// base class to make his/her own digi class. Two virtual method Draw()
// and Print() can be implemented if the user wants these functionarities.
//  If a concrete digi class is used as a transient class, G4Allocator
// must be used.

class G4VDigi 
{

  public:
      G4VDigi();
      virtual ~G4VDigi();

      G4int operator==(const G4VDigi &right) const;

      virtual void Draw();
      virtual void Print();

};

#endif

