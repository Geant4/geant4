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
// $Id: RunAction.hh,v 1.5 2004/06/21 10:57:11 maire Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;
class DetectorConstruction;
class PrimaryGeneratorAction;
class HistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:
    RunAction(DetectorConstruction*, PrimaryGeneratorAction*, HistoManager*);
   ~RunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);

    void AddEnergy (G4double edep)
                   {EnergyDeposit += edep; EnergyDeposit2 += edep*edep;};

    void AddTrakLenCharg (G4double length)
                 {TrakLenCharged += length; TrakLenCharged2 += length*length;};

    void AddTrakLenNeutr (G4double length)
                 {TrakLenNeutral += length; TrakLenNeutral2 += length*length;};

    void AddMscProjTheta (G4double theta)
                 {MscProjecTheta += theta;  MscProjecTheta2 += theta*theta;};

    void CountStepsCharg (G4int nSteps)
                 {nbStepsCharged += nSteps; nbStepsCharged2 += nSteps*nSteps;};

    void CountStepsNeutr (G4int nSteps)
                 {nbStepsNeutral += nSteps; nbStepsNeutral2 += nSteps*nSteps;};

    void CountParticles (G4ParticleDefinition* part)
                 {     if (part == G4Gamma::Gamma())       nbGamma++ ;
		  else if (part == G4Electron::Electron()) nbElect++ ;
		  else if (part == G4Positron::Positron()) nbPosit++ ; };

    void CountTransmit (G4int flag)
                 {     if (flag == 1)  Transmit[0]++;
		  else if (flag == 2) {Transmit[0]++; Transmit[1]++; }};

    void CountReflect (G4int flag)
                 {     if (flag == 1)  Reflect[0]++;
		  else if (flag == 2) {Reflect[0]++; Reflect[1]++; }};

    G4double ComputeMscHighland();

  private:
    G4double EnergyDeposit,  EnergyDeposit2;
    G4double TrakLenCharged, TrakLenCharged2;
    G4double TrakLenNeutral, TrakLenNeutral2;
    G4double nbStepsCharged, nbStepsCharged2;
    G4double nbStepsNeutral, nbStepsNeutral2;
    G4double MscProjecTheta, MscProjecTheta2;
    G4int    nbGamma, nbElect, nbPosit;
    G4int    Transmit[2],   Reflect[2];

    DetectorConstruction*   detector;
    PrimaryGeneratorAction* primary;
    HistoManager*           histoManager;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

