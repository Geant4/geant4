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
// -------------------------------------------------------------------
// $Id$
// -------------------------------------------------------------------

#ifndef RunAction_h
#define RunAction_h 1

#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>

#include "globals.hh"
#include "G4UserRunAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;

class RunAction : public G4UserRunAction
{
public:
  
  RunAction(DetectorConstruction*, PrimaryGeneratorAction*, HistoManager*);
  ~RunAction();

  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);
    
  G4int GetRow(){return row;}
  void AddRow() {row=row+1;}
  
  void AddToXVector(float v) {xVector(row)=v;}
  void AddToYVector(float v) {yVector(row)=v;}
  void AddToThetaVector(float v) {thetaVector(row)=v;}
  void AddToPhiVector(float v) {phiVector(row)=v;}
    
private:

  DetectorConstruction* detector;    
  PrimaryGeneratorAction* primary;   
  HistoManager* hist;
    
  // Matrix handling 
  
  G4int row;
  
  CLHEP::HepVector xVector;
  CLHEP::HepVector yVector;
  CLHEP::HepVector thetaVector;
  CLHEP::HepVector phiVector;
  CLHEP::HepMatrix beamMatrix;
};

#endif
