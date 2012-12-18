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
// $Id: G4EvManMessenger.hh,v 1.6 2006-11-03 03:11:13 asaim Exp $
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
    G4UIcmdWithoutParameter* storeEvtCmd;
};

#endif

