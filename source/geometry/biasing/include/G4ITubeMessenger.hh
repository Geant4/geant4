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
// $Id: G4ITubeMessenger.hh,v 1.2 2002-07-12 10:40:43 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ITubeMessenger
//
// Class description:
//
// a messenger createing  G4ISingleTubeMessenger which allow to
// create comands for a named cell: 
// /imp/cell/create <cellname>
// The G4ITubeFactory is used by G4ISingleTubeMessenger to create the 
// cells.
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4ITubeMessenger_hh
#define G4ITubeMessenger_hh G4ITubeMessenger_hh

#include "G4UImessenger.hh"
#include "g4std/map"

class G4UIcmdWithAString;
class G4ISingleTubeMessenger;
class G4ITubeFactory;

typedef G4std::map<G4String , G4ISingleTubeMessenger *> G4MapNameTubeMess;

class G4ITubeMessenger : public G4UImessenger {
public:
  G4ITubeMessenger(G4ITubeFactory *);
    // constructed with the  G4ITubeFactory to be messeged

  void SetNewValue(G4UIcommand * command, G4String newValue);

private:
  void Error(const G4String &m) {
    G4Exception("Error: G4ITubeMessenger" + m);
  }

  G4ITubeFactory *fITubeFactory;

  G4UIcmdWithAString *fCellCreateCmd;
  G4MapNameTubeMess fMapNameTubeMess;

};

#endif
