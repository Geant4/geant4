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
/// \file exoticphysics/monopole/include/RunAction.hh
/// \brief Definition of the RunAction class
//
// $Id: RunAction.hh 68036 2013-03-13 14:13:45Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class DetectorConstruction;
class RunActionMessenger;
class PrimaryGeneratorAction;
class G4Run;
class Histo;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
public:

  RunAction(DetectorConstruction*, PrimaryGeneratorAction*);
  virtual ~RunAction();

  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);

  void FillHisto(G4int id, G4double x, G4double weight = 1.0);
           
  //  G4double GetBinLength() {return binLength;};
  inline void SetBinSize(G4double size) { fBinLength =  size; }
  inline G4double GetOffsetX()          { return fOffsetX;} 

  inline void SetVerbose(G4int verbose) { fVerboseLevel = verbose;}
  inline G4int GetVerbose()             { return fVerboseLevel;}
    
  inline void AddProjRange (G4double x) { fProjRange += x; fProjRange2 += x*x; };
                   
private:  

  Histo*                  fHisto;    
  DetectorConstruction*   fDetector;
  PrimaryGeneratorAction* fKinematic;
  RunActionMessenger*     fRunActionMessenger;

  G4int                   fVerboseLevel;

  G4double                fBinLength;
  G4double                fOffsetX;
  G4double                fProjRange; 
  G4double                fProjRange2;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

