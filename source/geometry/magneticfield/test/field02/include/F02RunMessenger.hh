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
// $Id: F02RunMessenger.hh,v 1.3 2006-06-29 18:27:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F02RunMessenger_h
#define F02RunMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class F02RunAction;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class F02RunMessenger: public G4UImessenger
{
  public:

   F02RunMessenger(F02RunAction* );
  ~F02RunMessenger();

   void SetNewValue(G4UIcommand* ,G4String );

  private:

   F02RunAction*              runAction;
   
   G4UIdirectory*             plotDir;   

   G4UIcmdWithAString*        sethistNameCmd;

   G4UIcmdWithAnInteger*      setnbinStepCmd; 
   G4UIcmdWithADouble*        setSteplowCmd; 
   G4UIcmdWithADouble*        setStephighCmd; 

   G4UIcmdWithAnInteger*      setnbinEnCmd; 
   G4UIcmdWithADoubleAndUnit* setEnlowCmd; 
   G4UIcmdWithADoubleAndUnit* setEnhighCmd; 

   G4UIcmdWithAnInteger*      setnbinGammaCmd; 
   G4UIcmdWithADoubleAndUnit* setElowGammaCmd; 
   G4UIcmdWithADoubleAndUnit* setEhighGammaCmd; 

   G4UIcmdWithAnInteger*      setnbinTtCmd; 
   G4UIcmdWithADoubleAndUnit* setTtlowCmd; 
   G4UIcmdWithADoubleAndUnit* setTthighCmd; 

   G4UIcmdWithAnInteger*      setnbinTbCmd; 
   G4UIcmdWithADoubleAndUnit* setTblowCmd; 
   G4UIcmdWithADoubleAndUnit* setTbhighCmd; 

   G4UIcmdWithAnInteger*      setnbinTsecCmd; 
   G4UIcmdWithADoubleAndUnit* setTseclowCmd; 
   G4UIcmdWithADoubleAndUnit* setTsechighCmd; 

   G4UIcmdWithAnInteger*      setnbinTglobCmd; 
   G4UIcmdWithADoubleAndUnit* setTgloblowCmd; 
   G4UIcmdWithADoubleAndUnit* setTglobhighCmd; 

   G4UIcmdWithAnInteger*      setnbinRCmd; 
   G4UIcmdWithADoubleAndUnit* setRlowCmd; 
   G4UIcmdWithADoubleAndUnit* setRhighCmd; 

   G4UIcmdWithAnInteger*      setnbinThCmd; 
   G4UIcmdWithADoubleAndUnit* setThlowCmd; 
   G4UIcmdWithADoubleAndUnit* setThhighCmd; 

   G4UIcmdWithAnInteger*      setnbinThbackCmd; 
   G4UIcmdWithADoubleAndUnit* setThlowbackCmd; 
   G4UIcmdWithADoubleAndUnit* setThhighbackCmd; 

   G4UIcmdWithAnInteger*      setnbinzvertexCmd; 
   G4UIcmdWithADoubleAndUnit* setzlowCmd; 
   G4UIcmdWithADoubleAndUnit* setzhighCmd;
   
   G4UIdirectory*             RndmDir;
   G4UIcmdWithAnInteger*      RndmSaveCmd;    
   G4UIcmdWithAString*        RndmReadCmd;    
 
};

#endif

