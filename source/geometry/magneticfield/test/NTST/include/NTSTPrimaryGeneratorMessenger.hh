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
// $Id: NTSTPrimaryGeneratorMessenger.hh,v 1.2 2003-12-09 15:35:23 gunter Exp $
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








