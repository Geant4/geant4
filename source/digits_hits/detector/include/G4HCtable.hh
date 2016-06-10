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
// $Id: G4HCtable.hh 67992 2013-03-13 10:59:57Z gcosmo $
//

#ifndef G4HCtable_H
#define G4HCtable_H 1

class G4VSensitiveDetector;
#include "globals.hh"
//#include "g4rw/tvordvec.h"
#include <vector>

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
    G4int GetCollectionID(G4String HCname) const;
    G4int GetCollectionID(G4VSensitiveDetector* aSD) const;

  private:
    std::vector<G4String> SDlist;
    std::vector<G4String> HClist;

  public:
    inline G4int entries() const
    { return HClist.size(); }
    inline G4String GetSDname(G4int i) const
    {
      if(i<0||i>entries()) return "***Not Defined***";
      return SDlist[i];
    }
    inline G4String GetHCname(G4int i) const
    {
      if(i<0||i>entries()) return "***Not Defined***";
      return HClist[i];
    }

};

#endif

