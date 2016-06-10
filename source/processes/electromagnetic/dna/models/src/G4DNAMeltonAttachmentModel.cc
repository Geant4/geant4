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
// $Id: G4DNAMeltonAttachmentModel.cc 85244 2014-10-27 08:24:13Z gcosmo $
//

// Created by Z. Francis

#include "G4DNAMeltonAttachmentModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAMeltonAttachmentModel::G4DNAMeltonAttachmentModel(const G4ParticleDefinition*,
                                                       const G4String& nam) :
    G4VEmModel(nam), isInitialised(false)
{
//    nistwater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
  fpWaterDensity = 0;

  lowEnergyLimit = 4 * eV;
  lowEnergyLimitOfModel = 4 * eV;
  highEnergyLimit = 13 * eV;
  SetLowEnergyLimit(lowEnergyLimit);
  SetHighEnergyLimit(highEnergyLimit);

  verboseLevel = 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if (verboseLevel > 0)
  {
    G4cout << "Melton Attachment model is constructed " << G4endl<< "Energy range: "
    << lowEnergyLimit / eV << " eV - "
    << highEnergyLimit / eV << " eV"
    << G4endl;
  }
  fParticleChangeForGamma = 0;
  fDissociationFlag = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAMeltonAttachmentModel::~G4DNAMeltonAttachmentModel()
{
  // For total cross section

  std::map<G4String, G4DNACrossSectionDataSet*, std::less<G4String> >::iterator pos;

  for (pos = tableData.begin(); pos != tableData.end(); ++pos)
  {
    G4DNACrossSectionDataSet* table = pos->second;
    delete table;
  }

  // For final state

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAMeltonAttachmentModel::Initialise(const G4ParticleDefinition* /*particle*/,
                                            const G4DataVector& /*cuts*/)
{

  if (verboseLevel > 3) G4cout
      << "Calling G4DNAMeltonAttachmentModel::Initialise()" << G4endl;

  // Energy limits

  if (LowEnergyLimit() < lowEnergyLimit)
  {
    G4cout << "G4DNAMeltonAttachmentModel: low energy limit increased from " <<
    LowEnergyLimit()/eV << " eV to " << lowEnergyLimit/eV << " eV" << G4endl;
    SetLowEnergyLimit(lowEnergyLimit);
  }

  if (HighEnergyLimit() > highEnergyLimit)
  {
    G4cout << "G4DNAMeltonAttachmentModel: high energy limit decreased from " <<
    HighEnergyLimit()/eV << " eV to " << highEnergyLimit/eV << " eV" << G4endl;
    SetHighEnergyLimit(highEnergyLimit);
  }

  // Reading of data files

  G4double scaleFactor = 1e-18*cm*cm;

  G4String fileElectron("dna/sigma_attachment_e_melton");

  G4ParticleDefinition* electronDef = G4Electron::ElectronDefinition();
  G4String electron;

  // ELECTRON

  // For total cross section

  electron = electronDef->GetParticleName();

  tableFile[electron] = fileElectron;

  G4DNACrossSectionDataSet* tableE =
      new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
  tableE->LoadData(fileElectron);
  tableData[electron] = tableE;

  //

  if (verboseLevel > 2)
  G4cout << "Loaded cross section data for Melton Attachment model" << G4endl;

  if( verboseLevel>0 )
  {
    G4cout << "Melton Attachment model is initialized " << G4endl
    << "Energy range: "
    << LowEnergyLimit() / eV << " eV - "
    << HighEnergyLimit() / eV << " eV"
    << G4endl;
  }
  // Initialize water density pointer
  fpWaterDensity = G4DNAMolecularMaterial::Instance()->
      GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));

  if (isInitialised)
  { return;}
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4DNAMeltonAttachmentModel::CrossSectionPerVolume(const G4Material* material,
                                                  const G4ParticleDefinition* particleDefinition,
                                                  G4double ekin,
                                                  G4double,
                                                  G4double)
{
  if (verboseLevel > 3) G4cout
      << "Calling CrossSectionPerVolume() of G4DNAMeltonAttachmentModel"
      << G4endl;

  // Calculate total cross section for model

  G4double sigma=0;

  G4double waterDensity = (*fpWaterDensity)[material->GetIndex()];

  if(waterDensity!= 0.0)
  //  if (material == nistwater || material->GetBaseMaterial() == nistwater)
  {
    const G4String& particleName = particleDefinition->GetParticleName();

    if (ekin >= lowEnergyLimit && ekin < highEnergyLimit)
    {

      std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
      pos = tableData.find(particleName);

      if (pos != tableData.end())
      {
        G4DNACrossSectionDataSet* table = pos->second;
        if (table != 0)
        {
          sigma = table->FindValue(ekin);
        }
      }
      else
      {
        G4Exception("G4DNAMeltonAttachmentModel::ComputeCrossSectionPerVolume",
                    "em0002",
            FatalException,"Model not applicable to particle type.");
      }
    }

    if (verboseLevel > 2)
    {
      G4cout << "__________________________________" << G4endl;
      G4cout << "=== G4DNAMeltonAttachmentModel - XS INFO START" << G4endl;
      G4cout << "--- Kinetic energy(eV)=" << ekin/eV
          << " particle : " << particleDefinition->GetParticleName() << G4endl;
      G4cout << "--- Cross section per water molecule (cm^2)="
          << sigma/cm/cm << G4endl;
      G4cout << "--- Cross section per water molecule (cm^-1)="
          << sigma*waterDensity/(1./cm) << G4endl;
      // G4cout << "--- Cross section per water molecule (cm^-1)="
      //        << sigma*material->GetAtomicNumDensityVector()[1]/(1./cm)
      //        << G4endl;
      G4cout << "--- G4DNAMeltonAttachmentModel - XS INFO END" << G4endl;
    }

  } // if water

  return sigma*waterDensity;
//    return sigma*material->GetAtomicNumDensityVector()[1];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAMeltonAttachmentModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
                                                   const G4MaterialCutsCouple* /*couple*/,
                                                   const G4DynamicParticle* aDynamicElectron,
                                                   G4double,
                                                   G4double)
{

  if (verboseLevel > 3) G4cout
      << "Calling SampleSecondaries() of G4DNAMeltonAttachmentModel" << G4endl;

      // Electron is killed

      G4double electronEnergy0 = aDynamicElectron->GetKineticEnergy();
      fParticleChangeForGamma->SetProposedKineticEnergy(0.);
      fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(electronEnergy0);

      if(fDissociationFlag)
      {
        G4DNAChemistryManager::Instance()->CreateWaterMolecule(eDissociativeAttachment,-1,
            fParticleChangeForGamma->GetCurrentTrack());
      }
      return;
    }
