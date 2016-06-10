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
// $Id: G4MatScanMessenger.hh 66892 2013-01-17 10:57:59Z gunter $
//
//

// class description:
//


#ifndef G4MatScanMessenger_HH
#define G4MatScanMessenger_HH 1

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
    
    virtual G4String GetCurrentValue(G4UIcommand * command);
    virtual void SetNewValue(G4UIcommand * command,G4String newValue);

  private:
    G4MaterialScanner* theScanner;
    
    G4UIdirectory* msDirectory;
    G4UIcmdWithoutParameter* scanCmd;
    G4UIcommand* thetaCmd;
    G4UIcommand* phiCmd;
    G4UIcommand* singleCmd;
    G4UIcmdWith3Vector* single2Cmd;
    G4UIcmdWithABool* regSenseCmd;
    G4UIcmdWithAString* regionCmd;
    G4UIcmdWith3VectorAndUnit* eyePosCmd;
};

#endif



