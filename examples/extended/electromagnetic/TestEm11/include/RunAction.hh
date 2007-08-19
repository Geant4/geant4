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
// $Id: RunAction.hh,v 1.3 2007-08-19 20:52:53 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;
class DetectorConstruction;
class PhysicsList;
class PrimaryGeneratorAction;
class HistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:
    RunAction(DetectorConstruction*, PhysicsList*,PrimaryGeneratorAction*,
              HistoManager*);
   ~RunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);
    
    void AddEdep (G4double e)        { Edeposit  += e; Edeposit2  += e*e;};
    void AddTrackLength (G4double t) { trackLen  += t; trackLen2  += t*t;};
    void AddProjRange   (G4double x) { projRange += x; projRange2 += x*x;};
    void AddStepSize    (G4int nb, G4double s)
                                     { nbOfSteps += nb; nbOfSteps2 += nb*nb;
                                       stepSize  += s ; stepSize2  += s*s;  };
    void AddTrackStatus (G4int i)    { status[i]++ ;};			           
    
  private:
    DetectorConstruction*   detector;
    PhysicsList*            physics;
    PrimaryGeneratorAction* kinematic;
    HistoManager*           histoManager;
    
    G4double                Edeposit,  Edeposit2;
    G4double                trackLen,  trackLen2;
    G4double                projRange, projRange2;
    G4int                   nbOfSteps, nbOfSteps2;
    G4double                stepSize,  stepSize2;
    G4int                   status[3];
    
    G4double                csdaRange;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

