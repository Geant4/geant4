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
// $Id: G4DCtable.hh,v 1.6 2001-07-11 10:08:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4DCtable_H
#define G4DCtable_H 1

#include "globals.hh"
//#include "g4rw/tvordvec.h"
#include "g4std/vector"

// class description:
//
//  This class is used by G4DigiManager for book keeping the
// digitizer modules and digits collections. The order of
// digi collections stored in G4DCofThisEvent is same as the
// order of DClist. 
//  The order may vary from run to run, if the user adds/changes
// some of his/her digitizer modules.
//  In case user wants to make G4Run object persistent, this
// G4DCtable class object should be copied and stored with
// G4Run object.

class G4DCtable
{
  public:
    G4DCtable();
    ~G4DCtable();

  public:
    G4int Registor(G4String SDname,G4String DCname);
    G4int GetCollectionID(G4String DCname);

  private:
    G4std::vector<G4String> DMlist;
    G4std::vector<G4String> DClist;

  public:
    inline G4int entries() const
    { return DClist.size(); }

};

#endif

