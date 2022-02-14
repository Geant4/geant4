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
// G4MatScanMessenger
//
// Class description:
//
// Messenger class for materials scanner.

// Author: M.Asai, May 2006
// --------------------------------------------------------------------
#ifndef G4MatScanMessenger_hh
#define G4MatScanMessenger_hh 1

#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcommand;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4MaterialScanner;

class G4MatScanMessenger : public G4UImessenger
{
  public:

    G4MatScanMessenger(G4MaterialScanner* p1);
    virtual ~G4MatScanMessenger();

    virtual G4String GetCurrentValue(G4UIcommand* command);
    virtual void SetNewValue(G4UIcommand* command, G4String newValue);

  private:

    G4MaterialScanner* theScanner = nullptr;

    G4UIdirectory* msDirectory = nullptr;
    G4UIcmdWithoutParameter* scanCmd = nullptr;
    G4UIcommand* thetaCmd = nullptr;
    G4UIcommand* phiCmd = nullptr;
    G4UIcommand* singleCmd = nullptr;
    G4UIcmdWith3Vector* single2Cmd = nullptr;
    G4UIcmdWithABool* regSenseCmd = nullptr;
    G4UIcmdWithAString* regionCmd = nullptr;
    G4UIcmdWith3VectorAndUnit* eyePosCmd = nullptr;
};

#endif
