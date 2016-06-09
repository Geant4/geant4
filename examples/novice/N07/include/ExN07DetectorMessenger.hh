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
// $Id: ExN07DetectorMessenger.hh,v 1.3 2004/11/17 03:07:30 kurasige Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//

#ifndef ExN07DetectorMessenger_h
#define ExN07DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class ExN07DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;

class ExN07DetectorMessenger: public G4UImessenger
{
  public:
    ExN07DetectorMessenger(ExN07DetectorConstruction* );
    virtual ~ExN07DetectorMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    virtual G4String GetCurrentValue(G4UIcommand * command);
    
  private:
    ExN07DetectorConstruction* ExN07Detector;
    
    G4UIdirectory*             N07Dir;
    G4UIcmdWithAString*        AddMaterCmd;
    G4UIcmdWithAString*        AbsMaterCmd;
    G4UIcmdWithAString*        GapMaterCmd;
    G4UIcmdWithAnInteger*      numLayerCmd;
    G4UIcmdWithABool*          SerialCmd;
    G4UIcmdWithAnInteger*      verboseCmd;
};


#endif

