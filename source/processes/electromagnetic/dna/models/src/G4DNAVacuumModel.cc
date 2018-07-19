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
// $Id$
//

#include "G4DNAVacuumModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"

G4DNAVacuumModel::G4DNAVacuumModel(const G4String& applyToMaterial, const G4ParticleDefinition*,
                                   const G4String& nam)
    : G4VDNAModel(nam, applyToMaterial)
{
    verboseLevel = 0;

    if( verboseLevel>0 )
    {
        G4cout << "G4DNAVacuumModel is constructed " << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAVacuumModel::~G4DNAVacuumModel()
{ 
    if (verboseLevel > 3)
        G4cout << "Calling G4DNAVacuumModel::Initialise()" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAVacuumModel::Initialise(const G4ParticleDefinition* particle,
                                  const G4DataVector& /*cuts*/, G4ParticleChangeForGamma*)
{

    if (verboseLevel > 3)
        G4cout << "Calling G4DNAVacuumModel::Initialise()" << G4endl;

    EnableForMaterialAndParticle("G4_Galactic", particle->GetParticleName() );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAVacuumModel::CrossSectionPerVolume(const G4Material* /*material*/,
                                                 const G4String& /*materialName*/,
                                                 const G4ParticleDefinition* /*particleDefinition*/,
                                                 G4double /*ekin*/,
                                                 G4double /*emin*/,
                                                 G4double /*emax*/)
{
    if (verboseLevel > 3)
        G4cout << "Calling CrossSectionPerVolume() of G4DNAVacuumModel" << G4endl;

    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAVacuumModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
                                         const G4MaterialCutsCouple* /*couple*/,
                                         const G4String& /*materialName*/,
                                         const G4DynamicParticle* /*aDynamicParticle*/,
                                         G4ParticleChangeForGamma* /*particleChangeForGamma*/,
                                         G4double /*tmin*/,
                                         G4double /*tmax*/)
{

    if (verboseLevel > 3)
        G4cout << "Calling SampleSecondaries() of G4DNAVacuumModel" << G4endl;

}

