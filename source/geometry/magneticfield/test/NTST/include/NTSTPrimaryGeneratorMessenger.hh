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
// $Id: NTSTPrimaryGeneratorMessenger.hh,v 1.3 2006-06-29 18:25:51 gunter Exp $
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef NTSTPrimaryGeneratorMessenger_h
#define NTSTPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "G4VPrimaryGenerator.hh"
#include "globals.hh"

class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIdirectory;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class NTSTPrimaryGeneratorMessenger: public G4UImessenger
{
public:
  NTSTPrimaryGeneratorMessenger();
  ~NTSTPrimaryGeneratorMessenger();

  inline G4VPrimaryGenerator* GetGenerator(){return generator;}

  void SetNewValue(G4UIcommand*, G4String);
  inline const G4String& GetName() const {return Name;}
  inline const G4String& GetNames()const {return Names;}
  inline G4bool PrintState() const { return print; }
  inline void EnablePrinting() { print=true;}
  inline void DisablePrinting(){ print=false;}

private:
  G4UIdirectory*       GenDir;
  G4UIcmdWithAString*  ChooseCmd;
  G4UIcmdWithABool*    PrintCmd;
  G4String             Name;
  enum Choices         {GUN, EVT};
  G4String             Names;
  Choices              Choice;
  G4VPrimaryGenerator* generator;
  G4bool print;

};

#endif








