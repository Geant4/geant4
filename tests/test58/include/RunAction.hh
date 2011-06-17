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
// RunAction.hh
// 
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
class PhysicsList;
class PrimaryGeneratorAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
public:
  RunAction(DetectorConstruction*,PhysicsList* ,PrimaryGeneratorAction*) ;
  ~RunAction();

public:
  void BeginOfRunAction(const G4Run*);
  void   EndOfRunAction(const G4Run*);

  void AddEnergy (G4double edep)
       { EnergyDeposit += edep; EnergyDeposit2 += edep*edep; };

  void AddNonIonEnergy (G4double enondep)
       { NonIonEnergyDeposit += enondep; 
         NonIonEnergyDeposit2 += enondep*enondep; };


   // sum of secondary Kinetic energy*L(T):
  void SumTL(G4double energyL){sum_TL +=energyL;  sum_TL2 +=energyL*energyL;};

  // sum of secondary Kinetic energy:
  void SumT(G4double energy){sum_T +=energy;  sum_T2 +=energy*energy; };

  //number of recoils
  void NumberRec(G4int i){ N_rec += i; };

  //number of steps
  void AddNumberOfSteps(G4int i ){ Nsteps+=double(i); 
				   Nsteps2+=double(i)*double(i);};

  //scattering angle	
  void AddTheta(G4double tet){ theta+=tet;
			       theta2+=tet*tet;};	
  // track length	
  void AddTrakLenPrim (G4double length)
          {TrakLenPrim += length; TrakLenPrim2 += length*length;};

    // Get Threshold energy for displacement
  void  GiveThreshold(G4double energy) {Th=energy;};

private:

  G4double EnergyDeposit,  EnergyDeposit2;
  G4double NonIonEnergyDeposit,  NonIonEnergyDeposit2;
  G4double sum_TL, sum_TL2;
  G4double sum_T, sum_T2;
  G4double Th; 
  G4int N_rec; 
  G4double Nsteps, Nsteps2;
  G4double theta, theta2;
  G4double TrakLenPrim, TrakLenPrim2;

  DetectorConstruction*   detector;
  PhysicsList*        physics;
  PrimaryGeneratorAction* primary;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

