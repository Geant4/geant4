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
// $Id: RunAction.hh,v 1.9 2009-01-22 17:41:43 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
    virtual ~RunAction();

    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);

    void AddEnergy (G4double edep)
                 {EnergyDeposit += edep; EnergyDeposit2 += edep*edep;};

    void AddTrakLenCharg (G4double length)
                 {TrakLenCharged += length; TrakLenCharged2 += length*length;};

    void AddTrakLenNeutr (G4double length)
                 {TrakLenNeutral += length; TrakLenNeutral2 += length*length;};

    void AddMscProjTheta (G4double theta)
                 {if (std::abs(theta) <= MscThetaCentral) { MscEntryCentral++;
		    MscProjecTheta += theta;  MscProjecTheta2 += theta*theta;}
		 };

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
    
    void AddEnergyLeak (G4double eleak, G4int index)
               { EnergyLeak[index] += eleak; EnergyLeak2[index] += eleak*eleak;};

  private:
    G4double EnergyDeposit,  EnergyDeposit2;
    G4double TrakLenCharged, TrakLenCharged2;
    G4double TrakLenNeutral, TrakLenNeutral2;
    G4double nbStepsCharged, nbStepsCharged2;
    G4double nbStepsNeutral, nbStepsNeutral2;
    G4double MscProjecTheta, MscProjecTheta2;
    G4double MscThetaCentral;
    
    G4int    nbGamma, nbElect, nbPosit;
    G4int    Transmit[2],   Reflect[2];
    G4int    MscEntryCentral;
    
    G4double EnergyLeak[2],  EnergyLeak2[2];

    DetectorConstruction*   detector;
    PrimaryGeneratorAction* primary;
    HistoManager*           histoManager;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

