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
// $Id: G4SDmessenger.hh 67992 2013-03-13 10:59:57Z gcosmo $
//

#ifndef G4SDmessenger_h
#define G4SDmessenger_h 1

#include "G4UImessenger.hh"

class G4SDManager;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

// class description:
//
//  This is a cncrete class of G4UImessenger which handles the commands for
// G4SDManager. This class has the following commands:
//   /hits/
//   /hits/list
//   /hits/activate
//   /hits/inactivate
//   /hts/verbose
//

class G4SDmessenger: public G4UImessenger
{
  public:
    G4SDmessenger(G4SDManager * SDManager);
    ~G4SDmessenger();
    void SetNewValue(G4UIcommand * command,G4String newValues);
  
  private:
    G4SDManager * fSDMan;
    G4UIdirectory* hitsDir;
    G4UIcmdWithoutParameter* listCmd;
    G4UIcmdWithAString* activeCmd;
    G4UIcmdWithAString* inactiveCmd;
    G4UIcmdWithAnInteger* verboseCmd;
};




#endif

