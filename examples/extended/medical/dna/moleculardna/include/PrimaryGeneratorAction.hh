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
//
/// \file PrimaryGeneratorAction.hh
/// \brief Header file for Primary Generator for Molecular DNA simulation

#ifndef MOLECULAR_PRIMARY_GENERATORACTION_HH
#define MOLECULAR_PRIMARY_GENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "PrimaryGeneratorMessenger.hh"

class G4GeneralParticleSource;
class PrimaryGeneratorSourceGRASCSV;
class PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
    public:
        PrimaryGeneratorAction();

        ~PrimaryGeneratorAction() override;

        void GeneratePrimaries(G4Event* anEvent) override;

        // Generate events from Phase Space file
        void SetInputFileName(const G4String& filename)
        {
          fMyInputFileName = filename;
        }

    private:
        G4GeneralParticleSource* fParticleGun = nullptr;

        // parameters for input from file
        G4String fMyInputFileName;
        static PrimaryGeneratorSourceGRASCSV* fPrimarySource;
        G4bool fFirstEvent;

        // class messenger
        PrimaryGeneratorMessenger*  fGunMessenger = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif  // MOLECULAR_PRIMARY_GENERATORACTION_HH
