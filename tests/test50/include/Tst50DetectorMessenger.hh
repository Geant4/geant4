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
// $Id: Tst50DetectorMessenger.hh,v 1.8 2006-06-29 22:04:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// author: Susanna Guatelli (guatelli@ge.infn.it)
// 
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
//
// -------------------------------------------------------------------

#ifndef Tst50DetectorMessenger_h
#define Tst50DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;
class Tst50DetectorConstruction;

class Tst50DetectorMessenger: public G4UImessenger
{
public:
  Tst50DetectorMessenger(Tst50DetectorConstruction* );
  ~Tst50DetectorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  Tst50DetectorConstruction* detector;
    
  G4UIdirectory*             tst50Dir;
  G4UIcmdWithAString*        targetMaterialCmd;
  G4UIcmdWithADoubleAndUnit* targetThicknessCmd;
  G4UIcmdWithADoubleAndUnit* targetXDimensionCmd;
  G4UIcmdWithADoubleAndUnit* targetYDimensionCmd;             
  G4UIcmdWithABool* UseUserLimitCmd; 
  G4UIcmdWithADoubleAndUnit* setStepMaxCmd;
  G4UIcmdWithoutParameter*   updateCmd;
};
#endif

