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
// $Id: G4DNATransformElectronModel.cc 64057 2012-10-30 15:04:49Z gcosmo $
//
#include "G4DNATransformElectronModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Electron.hh"
#include "G4NistManager.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"

G4DNATransformElectronModel::G4DNATransformElectronModel(const G4ParticleDefinition*,
                                                         const G4String& nam):
    G4VEmModel(nam),fIsInitialised(false)
{
    fVerboseLevel = 0 ;
    SetLowEnergyLimit(0.*eV);
    SetHighEnergyLimit(0.025*eV);
    fParticleChangeForGamma = 0;
  //  fNistWater  = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    fpWaterDensity = 0;
    fpWaterDensity = 0;
    fEpsilon = 0.0001*eV;
}

//______________________________________________________________________
G4DNATransformElectronModel::~G4DNATransformElectronModel()
{}

//______________________________________________________________________
void G4DNATransformElectronModel::Initialise(const G4ParticleDefinition* particleDefinition,
                                             const G4DataVector&)
{
#ifdef G4VERBOSE
    if (fVerboseLevel)
        G4cout << "Calling G4DNATransformElectronModel::Initialise()" << G4endl;
#endif

    if (particleDefinition != G4Electron::ElectronDefinition())
    {
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "Attempting to calculate cross section for wrong particle";
        G4Exception("G4DNATransformElectronModel::CrossSectionPerVolume","G4DNATransformElectronModel001",
                    FatalErrorInArgument,exceptionDescription);
        return;
    }

    // Initialize water density pointer
    fpWaterDensity = G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));

    if(!fIsInitialised)
    {
        fIsInitialised = true;
        fParticleChangeForGamma = GetParticleChangeForGamma();
    }
}

//______________________________________________________________________
G4double G4DNATransformElectronModel::CrossSectionPerVolume(const G4Material* material,
                                                            const G4ParticleDefinition*,
                                                            G4double ekin,
                                                            G4double,
                                                            G4double)
{
#if G4VERBOSE
    if (fVerboseLevel > 1)
        G4cout << "Calling CrossSectionPerVolume() of G4DNATransformElectronModel" << G4endl;
#endif

    if(ekin - fEpsilon > HighEnergyLimit())
    {
        return 0.0;
    }

    G4double waterDensity = (*fpWaterDensity)[material->GetIndex()];

    if(waterDensity!= 0.0)
   //  if (material == nistwater || material->GetBaseMaterial() == nistwater)
    {
        if (ekin - fEpsilon <= HighEnergyLimit())
        {
            return DBL_MAX;
        }
    }

    return 0.0 ;
}

//______________________________________________________________________
void G4DNATransformElectronModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
                                                    const G4MaterialCutsCouple*,
                                                    const G4DynamicParticle* particle,
                                                    G4double,
                                                    G4double)
{
#if G4VERBOSE
    if (fVerboseLevel)
        G4cout << "Calling SampleSecondaries() of G4DNATransformElectronModel" << G4endl;
#endif

    G4double k = particle->GetKineticEnergy();

//    if (k - fEpsilon <= HighEnergyLimit())
//    {
        const G4Track * track = fParticleChangeForGamma->GetCurrentTrack();
        G4DNAChemistryManager::Instance()->CreateSolvatedElectron(track);
        fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
        fParticleChangeForGamma->ProposeLocalEnergyDeposit(k);
//    }
}
