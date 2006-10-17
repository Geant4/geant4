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
//    **************************************
//    *                                    *
//    *     CellDetectorMessenger.hh       *
//    *                                    *
//    **************************************
//
// 
// author: Susanna Guatelli (guatelli@ge.infn.it)
//	   Barbara Mascialino (Barbara.Mascialino@ge.infn.it)
// 
// History:
// -----------
// 20 September 2006   S. Guatelli, B. Mascialino   1st implementation
//
// -------------------------------------------------------------------

#ifndef CellDetectorMessenger_h
#define CellDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;
class CellDetectorConstruction;

class CellDetectorMessenger: public G4UImessenger
{
public:
  CellDetectorMessenger(CellDetectorConstruction* );
  ~CellDetectorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  CellDetectorConstruction* detector;
    
  G4UIdirectory*             cellDir;
  G4UIcmdWithAString*        targetMaterialCmd;
  G4UIcmdWithADoubleAndUnit* targetThicknessCmd;
  G4UIcmdWithADoubleAndUnit* targetXDimensionCmd;
  G4UIcmdWithADoubleAndUnit* targetYDimensionCmd;             
  G4UIcmdWithABool* UseUserLimitCmd; 
  G4UIcmdWithADoubleAndUnit* setStepMaxCmd;
  G4UIcmdWithoutParameter*   updateCmd;
};
#endif

