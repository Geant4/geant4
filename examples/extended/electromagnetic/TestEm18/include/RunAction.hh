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
// $Id: RunAction.hh,v 1.1 2007-02-13 17:57:20 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
    RunAction(DetectorConstruction*, PrimaryGeneratorAction*, HistoManager*);
   ~RunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);

    void AddEnergyDeposit (G4double edep)
                   {energyDeposit += edep;};

    void AddTrackLength (G4double step)
                 {trackLength += step; nbSteps++;};
		 
    void AddChargedSecondary (G4double ekin)
                 {energyCharged += ekin; nbCharged++;
		  if (ekin<emin[0]) emin[0] = ekin;
		  if (ekin>emax[0]) emax[0] = ekin;
		 };
		 
    void AddNeutralSecondary (G4double ekin)
                 {energyNeutral += ekin; nbNeutral++;
		  if (ekin<emin[1]) emin[1] = ekin;
		  if (ekin>emax[1]) emax[1] = ekin;
		 };
		
  public:
    G4double GetEnergyFromRestrictedRange
             (G4double,G4ParticleDefinition*,G4Material*,G4double);
	     	  
    G4double GetEnergyFromCSDARange
             (G4double,G4ParticleDefinition*,G4Material*,G4double);		 
		 
  private:
    G4double energyDeposit;
    G4double trackLength;
    G4double energyCharged, energyNeutral;
    G4double emin[2], emax[2];
    
    G4long   nbSteps;
    G4int    nbCharged, nbNeutral;

    DetectorConstruction*   detector;
    PrimaryGeneratorAction* primary;
    HistoManager*           histoManager;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

