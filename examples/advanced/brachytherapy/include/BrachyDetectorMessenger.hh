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
//    *****************************************
//    *                                       *
//    *      BrachyDetectrorMessenger.hh      *
//    *                                       *
//    *****************************************
//
// $Id: BrachyDetectorMessenger.hh,v 1.4 2002-11-18 15:18:36 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef BrachyDetectorMessenger_h
#define BrachyDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class BrachyDetectorConstruction;
class BrachyFactoryIr;
class BrachyRunAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class BrachyDetectorMessenger: public G4UImessenger
{
  public:
    BrachyDetectorMessenger(BrachyDetectorConstruction* );
   ~BrachyDetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
  
  private:

BrachyDetectorConstruction*  Detector;
G4UIcmdWithAString*        selDetCmd;
G4UIcmdWithAString*        switchCmd;
G4UIcmdWithAString*        cleanCmd;
G4UIdirectory*             detDir;
G4UIdirectory*           mydetDir;
G4UIcmdWithAString*        AbsMaterCmd;
G4int  detectorChoice; 
G4int flag;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

