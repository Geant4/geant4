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
// $Id: G4UIdirectory.hh 74098 2013-09-22 16:05:49Z gcosmo $
//
//

#ifndef G4UIdirectory_H
#define G4UIdirectory_H 1

#include "G4UIcommand.hh"

// class description:
//  A concrete class of G4UIcommand. This class defines a command
// directory which can have commands.
//  General information of G4UIcommand is given in G4UIcommand.hh.

class G4UIdirectory : public G4UIcommand
{
  public: // with description
    G4UIdirectory(char * theCommandPath,G4bool commandsToBeBroadcasted = true);
    G4UIdirectory(const char * theCommandPath,G4bool commandsToBeBroadcasted = true);
    // Constructors. The argument is a full path directory which
    // starts and ends with "/".
};

#endif
