// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelRunAction.hh,v 1.3 2000-12-06 16:53:13 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayTelRunAction  ------
//           by  R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef GammaRayTelRunAction_h
#define GammaRayTelRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

#ifdef G4ANALYSIS_USE
#include "GammaRayTelAnalysisManager.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;

class GammaRayTelRunAction : public G4UserRunAction
{
public:
#ifdef G4ANALYSIS_USE
  GammaRayTelRunAction(GammaRayTelAnalysisManager* analysisMgr);
#else
  GammaRayTelRunAction();
#endif
  ~GammaRayTelRunAction();
  
public:
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);

private:
#ifdef G4ANALYSIS_USE
    GammaRayTelAnalysisManager* analysisManager;
#endif  
};

#endif



