// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20PrimaryGeneratorMessenger.hh,v 1.1 2001-05-24 19:49:21 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ Tst20PrimaryGeneratorMessenger  ------
// ************************************************************

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst20PrimaryGeneratorMessenger_h
#define Tst20PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class Tst20PrimaryGeneratorAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst20PrimaryGeneratorMessenger: public G4UImessenger
{
public:
  Tst20PrimaryGeneratorMessenger(Tst20PrimaryGeneratorAction*);
  ~Tst20PrimaryGeneratorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  Tst20PrimaryGeneratorAction* Tst20Action; 
  G4UIcmdWithAString*          RndmCmd;
  G4UIcmdWithAnInteger*        SourceTypeCmd;
  G4UIcmdWithADoubleAndUnit*   VertexRadiusCmd;
  G4UIcmdWithAnInteger*        SpectrumTypeCmd;
};

#endif


