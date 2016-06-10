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
// $Id: G4DMmessenger.hh 80987 2014-05-19 10:50:22Z gcosmo $
//

#ifndef G4DMmessenger_h
#define G4DMmessenger_h 1

#include "G4UImessenger.hh"

class G4DigiManager;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;

// class description:
//
//  This class is a concrete class of G4UImessenger which manages
// commands for G4DigiManager. Commands defined in this class are
//    /digi/
//    /digi/List
//    /digi/Digitize
//    /digi/Verbose
// These commands are available only if the user creates his/her
// digitizer module(s).
//

class G4DMmessenger: public G4UImessenger
{
  public:
    G4DMmessenger(G4DigiManager * DigiManager);
    ~G4DMmessenger();
    void SetNewValue(G4UIcommand * command,G4String newValues);
  
  private:
    G4DigiManager * fDMan;
    G4UIdirectory* digiDir;
    G4UIcmdWithoutParameter* listCmd;
    G4UIcmdWithAString* digiCmd;
    G4UIcmdWithAnInteger* verboseCmd;
};




#endif

