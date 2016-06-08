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
// $Id: G4SDmessenger.hh,v 1.3.2.1 2001/06/28 19:07:48 gunter Exp $
// GEANT4 tag $Name:  $
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

