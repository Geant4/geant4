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
//
// $Id: G4DCtable.hh 80987 2014-05-19 10:50:22Z gcosmo $
//

#ifndef G4DCtable_H
#define G4DCtable_H 1

class G4VDigitizerModule;
#include "globals.hh"
//#include "g4rw/tvordvec.h"
#include <vector>

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
    G4int GetCollectionID(G4String DCname) const;
    G4int GetCollectionID(G4VDigitizerModule* aDM) const;

  private:
    std::vector<G4String> DMlist;
    std::vector<G4String> DClist;

  public:
    inline G4int entries() const
    { return DClist.size(); }
    inline G4String GetDMname(G4int i) const
    {
      if(i<0||i>entries()) return "***Not Defined***";
      return DMlist[i];
    }
    inline G4String GetDCname(G4int i) const
    {
      if(i<0||i>entries()) return "***Not Defined***";
      return DClist[i];
    }

};

#endif

