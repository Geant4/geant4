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
// $Id: G4ITubeMessenger.cc,v 1.2 2002-07-12 10:40:44 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ITubeMessenger.cc
//
// ----------------------------------------------------------------------

#include "G4ITubeMessenger.hh"
#include "G4ITubeFactory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ISingleTubeMessenger.hh"


G4ITubeMessenger::G4ITubeMessenger(G4ITubeFactory *itfact) :
  fITubeFactory(itfact)
{  
  fCellCreateCmd = new G4UIcmdWithAString("/imp/cell/create", this);
}

void G4ITubeMessenger::
SetNewValue(G4UIcommand * command, G4String newValue){
  if (command==fCellCreateCmd) {
    G4MapNameTubeMess::iterator it = fMapNameTubeMess.find(newValue);
    if (it!=fMapNameTubeMess.end()) {
      Error("cell: " + newValue + ", already exists!");
    }
    fMapNameTubeMess[newValue] = new G4ISingleTubeMessenger(newValue, 
							fITubeFactory);
  }
}

