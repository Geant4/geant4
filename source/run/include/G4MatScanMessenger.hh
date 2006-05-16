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
// $Id: G4MatScanMessenger.hh,v 1.2 2006-05-16 21:57:14 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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



