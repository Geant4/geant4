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
// $Id: Em3DetectorMessenger.hh,v 1.5 2001-10-22 10:58:48 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Em3DetectorMessenger_h
#define Em3DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em3DetectorConstruction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Em3DetectorMessenger: public G4UImessenger
{
  public:
    Em3DetectorMessenger(Em3DetectorConstruction* );
   ~Em3DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Em3DetectorConstruction* Em3Detector;
    
    G4UIdirectory*             Em3detDir;

    G4UIcmdWithADoubleAndUnit* SizeYZCmd;
    G4UIcmdWithAnInteger*      NbLayersCmd;
    G4UIcmdWithAnInteger*      NbAbsorCmd;
    G4UIcommand*               AbsorCmd;        
    G4UIcmdWithADoubleAndUnit* MagFieldCmd;
    G4UIcmdWithADoubleAndUnit* MaxStepCmd;
    G4UIcmdWithoutParameter*   UpdateCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

