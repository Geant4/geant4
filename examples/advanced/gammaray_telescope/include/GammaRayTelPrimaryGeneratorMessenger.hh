// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelPrimaryGeneratorMessenger.hh,v 1.2 2000-11-15 20:27:39 flongo Exp $
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
//      ------------ GammaRayTelPrimaryGeneratorMessenger  ------
//           by G.Santin, F.Longo & R.Giannitrapani (13 nov 2000) 
//
// ************************************************************

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef GammaRayTelPrimaryGeneratorMessenger_h
#define GammaRayTelPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class GammaRayTelPrimaryGeneratorAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelPrimaryGeneratorMessenger: public G4UImessenger
{
public:
  GammaRayTelPrimaryGeneratorMessenger(GammaRayTelPrimaryGeneratorAction*);
  ~GammaRayTelPrimaryGeneratorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  GammaRayTelPrimaryGeneratorAction* GammaRayTelAction; 
  G4UIcmdWithAString*          RndmCmd;
  G4UIcmdWithAnInteger*        SourceTypeCmd;
  G4UIcmdWithADoubleAndUnit*   VertexRadiusCmd;
  G4UIcmdWithAnInteger*        SpectrumTypeCmd;
};

#endif


