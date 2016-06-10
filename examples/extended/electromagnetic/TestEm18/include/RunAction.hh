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
/// \file electromagnetic/TestEm18/include/RunAction.hh
/// \brief Definition of the RunAction class
//
// $Id: RunAction.hh 66241 2012-12-13 18:34:42Z gunter $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;
class G4ParticleDefinition;
class G4Material;

class DetectorConstruction;
class PrimaryGeneratorAction;
class HistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:
    RunAction(DetectorConstruction*, PrimaryGeneratorAction*);
   ~RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    void AddEnergyDeposit (G4double edep)
                   {fEnergyDeposit += edep;};

    void AddTrackLength (G4double step)
                 {fTrackLength += step; fNbSteps++;};
                 
    void AddChargedSecondary (G4double ekin)
                 {fEnergyCharged += ekin; fNbCharged++;
                  if (ekin<fEmin[0]) fEmin[0] = ekin;
                  if (ekin>fEmax[0]) fEmax[0] = ekin;
                 };
                 
    void AddNeutralSecondary (G4double ekin)
                 {fEnergyNeutral += ekin; fNbNeutral++;
                  if (ekin<fEmin[1]) fEmin[1] = ekin;
                  if (ekin>fEmax[1]) fEmax[1] = ekin;
                 };
                
  public:
    G4double GetEnergyFromRestrictedRange
             (G4double,G4ParticleDefinition*,G4Material*,G4double);
                       
    G4double GetEnergyFromCSDARange
             (G4double,G4ParticleDefinition*,G4Material*,G4double);                 
                 
  private:
    G4double fEnergyDeposit;
    G4double fTrackLength;
    G4double fEnergyCharged, fEnergyNeutral;
    G4double fEmin[2], fEmax[2];
    
    G4long   fNbSteps;
    G4int    fNbCharged, fNbNeutral;

    DetectorConstruction*   fDetector;
    PrimaryGeneratorAction* fPrimary;
    HistoManager*           fHistoManager;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

