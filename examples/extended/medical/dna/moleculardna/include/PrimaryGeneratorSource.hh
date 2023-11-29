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
/// \file PrimaryGeneratorSource.hh
/// \brief Header file for Primary Generator Source class for Molecular DNA simulation

#ifndef MOLECULAR_PRIMARY_GENERATORSOURCE_HH
#define MOLECULAR_PRIMARY_GENERATORSOURCE_HH

#include <list>
#include <fstream>

#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// - Abstract PrimaryGeneratorSource
//        - pure virtual functions
// - Implementations:
//        - PrimaryGeneratorSourceGRASCSV
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Class to handle Primary objects
class Primary
{
    public:
        // Primary class constructor
        Primary(){;}

        // Primary class destructor
        ~Primary() = default;

        void SetName(const G4String& name) { fName = name; }
        void SetName(G4int PDGnum)
        {
          if(PDGnum!=0)
          {
              G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
              G4ParticleDefinition* particle = pTable->FindParticle( PDGnum );
              fName = particle->GetParticleName();
          }
          else
          {
              fName="";
          }
        }
  
        G4String GetName() const { return fName; }

        void SetPosition(const G4ThreeVector& position) { fPosition = position; }
        G4ThreeVector GetPosition() const { return fPosition; }

        void SetMomentum(const G4ThreeVector& momentum) { fMomentum = momentum; }
        G4ThreeVector GetMomentum() const { return fMomentum; }

        void SetMomentumDirection(const G4ThreeVector& momentumDir) { fMomentumDir = momentumDir; }
        G4ThreeVector GetMomentumDirection() const { return fMomentumDir; }

        void SetEnergy(const G4double& energy) { fEnergy = energy; }
        G4double GetEnergy() const { return fEnergy; }

        G4ParticleDefinition* GetParticleDefinition()
        {
          G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
          G4ParticleDefinition* particle = pTable->FindParticle( fName.data() );
          return particle;
        }

        void Print()
        {
          G4cout << " *** " << fName << " *** " << G4endl;
          G4cout << " * Position: " << fPosition << G4endl;
          G4cout << " * Momentum Direction: " << fMomentumDir << G4endl;
          G4cout << " * Energy: " << fEnergy << G4endl;
          G4cout << " ***  *** ***" << G4endl;
        }

    private:

        G4String      fName;
        G4ThreeVector fPosition;
        G4ThreeVector fMomentum;
        G4ThreeVector fMomentumDir;
        G4double      fEnergy;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorSource
{
    public:

        // PrimaryGeneratorSource class destructor
        virtual ~PrimaryGeneratorSource() = default;
        
        // Define Primary class within PrimaryGeneratorSource
        virtual Primary* GetPrimary() = 0;
        virtual void SetBufferSize( G4int ) = 0;
        virtual void SetnParticles( G4int ) = 0;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

