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
// $Id: G4DMmessenger.hh,v 1.4 2001-07-11 10:08:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

