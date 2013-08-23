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
// Please cite the following paper if you use this software
// Nucl.Instrum.Meth.B260:20-27, 2007

#ifndef RunAction_h
#define RunAction_h 1

#include <CLHEP/Matrix/Vector.h>

#include "G4UserRunAction.hh"

#include "PrimaryGeneratorAction.hh"

class G4Run;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class RunAction : public G4UserRunAction
{
public:
  
  RunAction(DetectorConstruction*, PrimaryGeneratorAction*);
  ~RunAction();

  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);
    
  G4int GetRow(){return fRow;}
  void AddRow() {fRow=fRow+1;}
  
  void AddToXVector(float v) {fXVector(fRow)=v;}
  void AddToYVector(float v) {fYVector(fRow)=v;}
  void AddToThetaVector(float v) {fThetaVector(fRow)=v;}
  void AddToPhiVector(float v) {fPhiVector(fRow)=v;}
    
private:

  DetectorConstruction* fDetector;    
  PrimaryGeneratorAction* fPrimary;   
    
  // Matrix handling 
  
  G4int fRow;
  
  CLHEP::HepVector fXVector;
  CLHEP::HepVector fYVector;
  CLHEP::HepVector fThetaVector;
  CLHEP::HepVector fPhiVector;
  CLHEP::HepMatrix fBeamMatrix;
};

#endif
