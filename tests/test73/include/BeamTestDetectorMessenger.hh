//
//  BeamTestDetectorMessenger.hh
//  MSCTest
//
//  Created by Andrea Dotti on 3/7/12.
//  Copyright (c) 2012 CERN. All rights reserved.
//

#ifndef MSCTest_BeamTestDetectorMessenger_hh
#define MSCTest_BeamTestDetectorMessenger_hh

#include "globals.hh"
#include "G4UImessenger.hh"

class BeamTestDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;

class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;


class BeamTestDetectorMessenger: public G4UImessenger {
    
public:
    
    BeamTestDetectorMessenger( BeamTestDetectorConstruction* );
    ~BeamTestDetectorMessenger();
    
    void SetNewValue( G4UIcommand*, G4String );
    
private:
    
    BeamTestDetectorConstruction*    theDetector;
    G4UIdirectory*             theDetectorDir;
    G4UIcmdWithADoubleAndUnit* theChamberWidthCmd;
    G4UIcmdWithAnInteger*      theNumberOfChambersCmd;
    G4UIcmdWithADoubleAndUnit* theChamberSpacingCmd;    
    G4UIcmdWithoutParameter*   theUpdateCommand;
};



#endif
