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
// $Id: HadrontherapyDetectorMessenger.hh,v 1.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------

#ifndef HadrontherapyDetectorMessenger_h
#define HadrontherapyDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class HadrontherapyDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

// -------------------------------------------------------------------
class HadrontherapyDetectorMessenger: public G4UImessenger
{
  public:
    HadrontherapyDetectorMessenger(HadrontherapyDetectorConstruction* );
   ~HadrontherapyDetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
private:
  HadrontherapyDetectorConstruction* HadrontherapyDetector;
  
  G4UIdirectory*             HadronDir;
  G4UIdirectory*             detDir;
  G4UIcmdWithADoubleAndUnit* ModulatorAngleCmd;
};

// ----------------------------------------------------------------------

#endif

