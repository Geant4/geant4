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
/// \file PrimaryGeneratorSourceGRASCSV.hh
/// \brief Header file for Primary Generator source class implementation for GRAS CSV phase space for Molecular DNA simulation

#ifndef MOLECULAR_PRIMARY_GENERATORSOURCEGRASCSV_HH
#define MOLECULAR_PRIMARY_GENERATORSOURCEGRASCSV_HH

#include <list>
#include <fstream>

#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "PrimaryGeneratorSource.hh"

class PrimaryGeneratorSourceGRASCSV : public PrimaryGeneratorSource
{
    public:

        PrimaryGeneratorSourceGRASCSV(const G4String& filename);
        ~PrimaryGeneratorSourceGRASCSV();

        // Define Primary class within PrimaryGeneratorSourceGRASCSV
        Primary* GetPrimary();

        void SetBufferSize(G4int bufferSize){ fBufferSize = bufferSize; }
        void SetnParticles(G4int nParticles){ fnParticles = nParticles; }
        G4double RecomputeNParticles();

    private:
        std::ifstream       fInputFile;
        std::list<Primary*> fPrimaryList;
        G4int               fBufferSize = 100;
        G4int               fnParticles = 2;  // number of particles to be generated
        G4bool              fEndOfFile;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
