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
// Authors: S. Meylan and C. Villagrasa (IRSN, France)
// Models come from
// M. Bug et al, Rad. Phys and Chem. 130, 459-479 (2017)
//

#include "G4DNAPTBExcitationModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"

G4DNAPTBExcitationModel::G4DNAPTBExcitationModel(const G4String& applyToMaterial, const G4ParticleDefinition*,
                                                   const G4String& nam)
    : G4VDNAModel(nam, applyToMaterial)
{
    verboseLevel= 0;
    // Verbosity scale:
    // 0 = nothing
    // 1 = warning for energy non-conservation
    // 2 = details of energy budget
    // 3 = calculation of cross sections, file openings, sampling of atoms
    // 4 = entering in methods

    // initialisation of mean energy loss for each material
    tableMeanEnergyPTB["THF"] = 8.01*eV;
    tableMeanEnergyPTB["PY"] = 7.61*eV;
    tableMeanEnergyPTB["PU"] = 7.61*eV;
    tableMeanEnergyPTB["TMP"] = 8.01*eV;

    if( verboseLevel>0 )
    {
        G4cout << "PTB excitation model is constructed " << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAPTBExcitationModel::~G4DNAPTBExcitationModel()
{ 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAPTBExcitationModel::Initialise(const G4ParticleDefinition* particle,
                                         const G4DataVector& /*cuts*/, G4ParticleChangeForGamma*)
{
    if (verboseLevel > 3)
        G4cout << "Calling G4DNAPTBExcitationModel::Initialise()" << G4endl;

    G4double scaleFactor = 1e-16*cm*cm;
    G4double scaleFactorBorn = (1.e-22 / 3.343) * m*m;

    G4ParticleDefinition* electronDef = G4Electron::ElectronDefinition();

    //*******************************************************
    // Cross section data
    //*******************************************************

    if(particle == electronDef)
    {
        G4String particleName = particle->GetParticleName();

        AddCrossSectionData("THF",
                            particleName,
                            "dna/sigma_excitation_e-_PTB_THF",
                            scaleFactor);
        SetLowELimit("THF", particleName, 9.*eV);
        SetHighELimit("THF", particleName, 1.*keV);

        AddCrossSectionData("PY",
                            particleName,
                            "dna/sigma_excitation_e-_PTB_PY",
                            scaleFactor);
        SetLowELimit("PY", particleName, 9.*eV);
        SetHighELimit("PY", particleName, 1.*keV);

        AddCrossSectionData("PU",
                            particleName,
                            "dna/sigma_excitation_e-_PTB_PU",
                            scaleFactor);
        SetLowELimit("PU", particleName, 9.*eV);
        SetHighELimit("PU", particleName, 1.*keV);

        AddCrossSectionData("TMP",
                            particleName,
                            "dna/sigma_excitation_e-_PTB_TMP",
                            scaleFactor);
        SetLowELimit("TMP", particleName, 9.*eV);
        SetHighELimit("TMP", particleName, 1.*keV);

        AddCrossSectionData("G4_WATER",
                            particleName,
                            "dna/sigma_excitation_e_born",
                            scaleFactorBorn);
        SetLowELimit("G4_WATER", particleName, 9.*eV);
        SetHighELimit("G4_WATER", particleName, 1.*keV);

        // DNA materials
        //
        AddCrossSectionData("backbone_THF",
                            particleName,
                            "dna/sigma_excitation_e-_PTB_THF",
                            scaleFactor*33./30);
        SetLowELimit("backbone_THF", particleName, 9.*eV);
        SetHighELimit("backbone_THF", particleName, 1.*keV);

        AddCrossSectionData("cytosine_PY",
                            particleName,
                            "dna/sigma_excitation_e-_PTB_PY",
                            scaleFactor*42./30);
        SetLowELimit("cytosine_PY", particleName, 9.*eV);
        SetHighELimit("cytosine_PY", particleName, 1.*keV);

        AddCrossSectionData("thymine_PY",
                            particleName,
                            "dna/sigma_excitation_e-_PTB_PY",
                            scaleFactor*48./30);
        SetLowELimit("thymine_PY", particleName, 9.*eV);
        SetHighELimit("thymine_PY", particleName, 1.*keV);

        AddCrossSectionData("adenine_PU",
                            particleName,
                            "dna/sigma_excitation_e-_PTB_PU",
                            scaleFactor*50./44);
        SetLowELimit("adenine_PU", particleName, 9.*eV);
        SetHighELimit("adenine_PU", particleName, 1.*keV);

        AddCrossSectionData("guanine_PU",
                            particleName,
                            "dna/sigma_excitation_e-_PTB_PU",
                            scaleFactor*56./44);
        SetLowELimit("guanine_PU", particleName, 9.*eV);
        SetHighELimit("guanine_PU", particleName, 1.*keV);

        AddCrossSectionData("backbone_TMP",
                            particleName,
                            "dna/sigma_excitation_e-_PTB_TMP",
                            scaleFactor*33./50);
        SetLowELimit("backbone_TMP", particleName, 9.*eV);
        SetHighELimit("backbone_TMP", particleName, 1.*keV);
    }

    //*******************************************************
    // Load data
    //*******************************************************

    LoadCrossSectionData(particle->GetParticleName() );

    //*******************************************************
    // Verbose
    //*******************************************************

    if( verboseLevel>0 )
    {
        G4cout << "PTB excitation model is initialized " << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBExcitationModel::CrossSectionPerVolume(const G4Material* /*material*/,
                                                        const G4String& materialName,
                                                        const G4ParticleDefinition* particleDefinition,
                                                        G4double ekin,
                                                        G4double /*emin*/,
                                                        G4double /*emax*/)
{
    if (verboseLevel > 3)
        G4cout << "Calling CrossSectionPerVolume() of G4DNAPTBExcitationModel" << G4endl;

    // Get the name of the current particle
    G4String particleName = particleDefinition->GetParticleName();

    // initialise variables
    G4double lowLim = 0;
    G4double highLim = 0;
    G4double sigma=0;

    // Get the low energy limit for the current particle
    lowLim = GetLowELimit(materialName, particleName);

    // Get the high energy limit for the current particle
    highLim = GetHighELimit(materialName, particleName);

    // Check that we are in the correct energy range
    if (ekin >= lowLim && ekin < highLim)
    {
        // Get the map with all the data tables
        TableMapData* tableData = GetTableData();

        // Retrieve the cross section value
        sigma = (*tableData)[materialName][particleName]->FindValue(ekin);

        if (verboseLevel > 2)
        {
            G4cout << "__________________________________" << G4endl;
            G4cout << "°°° G4DNAPTBExcitationModel - XS INFO START" << G4endl;
            G4cout << "°°° Kinetic energy(eV)=" << ekin/eV << " particle : " << particleName << G4endl;
            G4cout << "°°° Cross section per "<< materialName <<" molecule (cm^2)=" << sigma/cm/cm << G4endl;
            G4cout << "°°° G4DNAPTBExcitationModel - XS INFO END" << G4endl;
        }

    }

    // Return the cross section value
    return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAPTBExcitationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
                                                const G4MaterialCutsCouple* /*couple*/,
                                                const G4String& materialName,
                                                const G4DynamicParticle* aDynamicParticle,
                                                G4ParticleChangeForGamma* particleChangeForGamma,
                                                G4double /*tmin*/,
                                                G4double /*tmax*/)
{
    if (verboseLevel > 3)
        G4cout << "Calling SampleSecondaries() of G4DNAPTBExcitationModel" << G4endl;

    // Get the incident particle kinetic energy
    G4double k = aDynamicParticle->GetKineticEnergy();

    if(materialName!="G4_WATER")
    {
        // Retrieve the excitation energy for the current material
        G4double excitationEnergy = tableMeanEnergyPTB[materialName];

        // Calculate the new energy of the particle
        G4double newEnergy = k - excitationEnergy;

        // Check that the new energy is above zero before applying it the particle.
        // Otherwise, do nothing.
        if (newEnergy > 0)
        {
            particleChangeForGamma->ProposeMomentumDirection(aDynamicParticle->GetMomentumDirection());
            particleChangeForGamma->SetProposedKineticEnergy(newEnergy);
            particleChangeForGamma->ProposeLocalEnergyDeposit(excitationEnergy);
        }
    }
    else
    {
        const G4String& particleName = aDynamicParticle->GetDefinition()->GetParticleName();

        G4int level = RandomSelectShell(k,particleName, materialName);
        G4double excitationEnergy = waterStructure.ExcitationEnergy(level);
        G4double newEnergy = k - excitationEnergy;

        if (newEnergy > 0)
        {
            particleChangeForGamma->ProposeMomentumDirection(aDynamicParticle->GetMomentumDirection());
            particleChangeForGamma->SetProposedKineticEnergy(newEnergy);
            particleChangeForGamma->ProposeLocalEnergyDeposit(excitationEnergy);
        }

        const G4Track * theIncomingTrack = particleChangeForGamma->GetCurrentTrack();
        G4DNAChemistryManager::Instance()->CreateWaterMolecule(eExcitedMolecule,
                                                               level,
                                                               theIncomingTrack);
    }
}
