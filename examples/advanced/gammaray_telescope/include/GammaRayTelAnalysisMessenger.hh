// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelAnalysisMessenger.hh,v 1.1 2000-12-06 16:53:13 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayTelAnalysysMessenger  ------
//           by R.Giannitrapani, F.Longo & G.Santin (03 dec 2000)
//
// ************************************************************
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
#ifdef G4ANALYSIS_USE
#ifndef GammaRayTelAnalysisMessenger_h
#define GammaRayTelAnalysisMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class GammaRayTelAnalysisManager;
class G4UIdirectory;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelAnalysisMessenger: public G4UImessenger
{
public:
  GammaRayTelAnalysisMessenger(GammaRayTelAnalysisManager* );
  ~GammaRayTelAnalysisMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  GammaRayTelAnalysisManager* GammaRayTelAnalysis;
  G4UIdirectory*              GammaRayTelAnalysisDir;
  
  G4UIcmdWithAString*        Histo1DDrawCmd;
  G4UIcmdWithAString*        Histo2DDrawCmd;
  G4UIcmdWithAString*        Histo1DSaveCmd;
  G4UIcmdWithAString*        Histo2DSaveCmd;
  G4UIcmdWithAString*        Histo2DModeCmd;
};
#endif
#endif










