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
// Contact authors: S. Meylan, C. Villagrasa
//
// email: sylvain.meylan@symalgo-tech.com, carmen.villagrasa@irsn.fr

#include "../include/G4DNADummyModel.hh"

#include "G4SystemOfUnits.hh"

G4DNADummyModel::G4DNADummyModel(const G4String& applyToMaterial, const G4ParticleDefinition* p, const G4String& nam, G4VEmModel* emModel)
    : G4VDNAModel(nam, applyToMaterial)
{
    fpEmModel = emModel;
    fpParticleDef = p;
}

G4DNADummyModel::~G4DNADummyModel()
{
    // There is no need to delete the model because it will be done in some G4 class.
    //if(fpEmModel) delete fpEmModel;
}

void G4DNADummyModel::Initialise(const G4ParticleDefinition* particle, const G4DataVector& v, G4ParticleChangeForGamma* changeForGamme)
{
    fMaterialMolPerVol = G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER") );

    fpEmModel->SetParticleChange(changeForGamme, nullptr);
    fpEmModel->Initialise(particle, v);

    // MatManagSys
    EnableForMaterialAndParticle("G4_WATER", fpParticleDef->GetParticleName() );
    SetLowELimit("G4_WATER", fpParticleDef->GetParticleName(), fpEmModel->LowEnergyLimit() );
    SetHighELimit("G4_WATER",fpParticleDef->GetParticleName(), fpEmModel->HighEnergyLimit() );
}

G4double G4DNADummyModel::CrossSectionPerVolume(const G4Material* material, const G4String& /*materialName*/, const G4ParticleDefinition* p, G4double ekin, G4double emin, G4double emax)
{
    G4double crossSectionTimesDensity = fpEmModel->CrossSectionPerVolume(material, p, ekin, emin, emax);
    G4double crossSection = crossSectionTimesDensity / GetNumMoleculePerVolumeUnitForMaterial(G4Material::GetMaterial("G4_WATER") );

    return crossSection;
}

void G4DNADummyModel::SampleSecondaries(std::vector<G4DynamicParticle*>* a, const G4MaterialCutsCouple* b, const G4String& /*materialName*/, const G4DynamicParticle* c, G4ParticleChangeForGamma* /*particleChangeForGamma*/, G4double tmin, G4double tmax)
{
    fpEmModel->SampleSecondaries(a, b, c, tmin, tmax);
}

G4double G4DNADummyModel::GetNumMoleculePerVolumeUnitForMaterial(const G4Material* mat)
{
    return fMaterialMolPerVol->at(mat->GetIndex() );
}


