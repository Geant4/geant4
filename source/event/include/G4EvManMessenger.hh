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
// $Id: G4EvManMessenger.hh,v 1.4 2001-07-11 09:58:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4EvManMessenger_h
#define G4EvManMessenger_h 1

#include "G4UImessenger.hh"
class G4EventManager;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;

// class description:
//
//  This is a concrete class of G4UImessenger which takes care of commands
// addressed to G4EventManager. Commands handled by this messenger are
//     /event/
//     /event/abort
//     /event/verbose
//

class G4EvManMessenger: public G4UImessenger
{
  public:
    G4EvManMessenger(G4EventManager * fEvMan);
    ~G4EvManMessenger();
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);
  private:
    G4EventManager * fEvManager;
    G4UIdirectory* eventDirectory;
    G4UIcmdWithoutParameter* abortCmd;
    G4UIcmdWithAnInteger* verboseCmd;
};

#endif

