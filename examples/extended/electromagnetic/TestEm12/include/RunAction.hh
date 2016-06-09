//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: RunAction.hh,v 1.1 2005/07/22 11:08:48 maire Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

