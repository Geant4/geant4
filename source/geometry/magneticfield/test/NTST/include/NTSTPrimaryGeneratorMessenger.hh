// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: NTSTPrimaryGeneratorMessenger.hh,v 1.1 2003-11-07 21:30:28 japost Exp $
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








