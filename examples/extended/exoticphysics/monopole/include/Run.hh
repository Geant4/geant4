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
/// \file exoticphysics/monopole/include/Run.hh
/// \brief Definition of the Run class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "globals.hh"
#include "g4root.hh"

class PrimaryGeneratorAction;
class DetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
public:
  Run(DetectorConstruction*, PrimaryGeneratorAction*);
  ~Run();
    
  virtual void Merge(const G4Run*);
  void EndOfRun(double binLength);   

  void FillHisto(G4int id, G4double x, G4double weight = 1.0);
           
  inline void SetVerbose(G4int verbose) { fVerboseLevel = verbose; }
  inline G4int GetVerbose() const       { return fVerboseLevel; }
  inline G4double GetOffsetX() const    { return fOffsetX; }

  inline void AddProjRange (G4double x) { fProjRange += x; 
                                          fProjRange2 += x*x; };

private:
  
  DetectorConstruction*   fDetector;
  PrimaryGeneratorAction* fPrimary;    
  G4AnalysisManager*      fAnalysisManager;

  G4int                   fNevt;

  G4double                fOffsetX;
  G4double                fProjRange; 
  G4double                fProjRange2;

  G4int                   fVerboseLevel;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

