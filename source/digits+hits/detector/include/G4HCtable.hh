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
// $Id: G4HCtable.hh,v 1.6 2001/07/11 09:58:42 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $
//

#ifndef G4HCtable_H
#define G4HCtable_H 1

#include "globals.hh"
//#include "g4rw/tvordvec.h"
#include "g4std/vector"

// class description:
//
//  This class is used by G4SDManager for book keeping the
// sensitive detector modules and hits collections. The order of
// hits collections stored in G4HCofThisEvent is same as the
// order of HClist. 
//  The order may vary from run to run, if the user adds/changes
// some of his/her sensitive detector modules.
//  In case user wants to make G4Run object persistent, this
// G4HCtable class object should be copied and stored with
// G4Run object.

class G4HCtable
{
  public:
    G4HCtable();
    ~G4HCtable();

  public:
    G4int Registor(G4String SDname,G4String HCname);
    G4int GetCollectionID(G4String HCname);

  private:
    G4std::vector<G4String> SDlist;
    G4std::vector<G4String> HClist;

  public:
    inline G4int entries() const
    { return HClist.size(); }

};

#endif

