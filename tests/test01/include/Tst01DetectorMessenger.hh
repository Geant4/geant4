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

#ifndef Tst01DetectorMessenger_h
#define Tst01DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst01DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

class Tst01DetectorMessenger: public G4UImessenger
{
  public:

    Tst01DetectorMessenger(Tst01DetectorConstruction* myDC) ;

    void SetNewValue(G4UIcommand* command, G4String newValues) ;

  private:

    Tst01DetectorConstruction* myDetector ;

    G4UIdirectory*      mydetDir ;

    G4UIcmdWithAString* selDetCmd ;
    G4UIcmdWithAString* switchCmd ;
    G4UIcmdWithAString* selMatCmd ;

  // Select CSG/Boolean

    G4UIcmdWithAString* selCSGcmd ;
    G4UIcmdWithAString* switchCSGcmd ;

  // Switch CSG/Boolean

    G4UIcmdWithAString* selBoolCmd ;
    G4UIcmdWithAString* switchBoolCmd ;

};

#endif

