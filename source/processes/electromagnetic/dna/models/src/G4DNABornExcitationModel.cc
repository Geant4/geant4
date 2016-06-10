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
// $Id: G4DNABornExcitationModel.cc 87137 2014-11-25 09:12:48Z gcosmo $
//

#include "G4DNABornExcitationModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNABornExcitationModel::G4DNABornExcitationModel(const G4ParticleDefinition*,
                                                   const G4String& nam) :
    G4VEmModel(nam), isInitialised(false)
{
  //  nistwater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
  fpMolWaterDensity = 0;

  verboseLevel = 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if (verboseLevel > 0)
  {
    G4cout << "Born excitation model is constructed " << G4endl;
  }
  fParticleChangeForGamma = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNABornExcitationModel::~G4DNABornExcitationModel()
{
  // Cross section

  std::map<G4String, G4DNACrossSectionDataSet*, std::less<G4String> >::iterator pos;
  for (pos = tableData.begin(); pos != tableData.end(); ++pos)
  {
    G4DNACrossSectionDataSet* table = pos->second;
    delete table;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNABornExcitationModel::Initialise(const G4ParticleDefinition* particle,
                                          const G4DataVector& /*cuts*/)
{

  if (verboseLevel > 3)
  G4cout << "Calling G4DNABornExcitationModel::Initialise()" << G4endl;

  G4String fileElectron("dna/sigma_excitation_e_born");
  G4String fileProton("dna/sigma_excitation_p_born");

  G4ParticleDefinition* electronDef = G4Electron::ElectronDefinition();
  G4ParticleDefinition* protonDef = G4Proton::ProtonDefinition();

  G4String electron;
  G4String proton;

  G4double scaleFactor = (1.e-22 / 3.343) * m*m;

  // *** ELECTRON

  electron = electronDef->GetParticleName();

  tableFile[electron] = fileElectron;

  lowEnergyLimit[electron] = 9. * eV;
  highEnergyLimit[electron] = 1. * MeV;

  // Cross section

  G4DNACrossSectionDataSet* tableE = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
  tableE->LoadData(fileElectron);

  tableData[electron] = tableE;

  // *** PROTON

  proton = protonDef->GetParticleName();

  tableFile[proton] = fileProton;

  lowEnergyLimit[proton] = 500. * keV;
  highEnergyLimit[proton] = 100. * MeV;

  // Cross section

  G4DNACrossSectionDataSet* tableP = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
  tableP->LoadData(fileProton);

  tableData[proton] = tableP;

  //

  if (particle==electronDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[electron]);
    SetHighEnergyLimit(highEnergyLimit[electron]);
  }

  if (particle==protonDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[proton]);
    SetHighEnergyLimit(highEnergyLimit[proton]);
  }

  if( verboseLevel>0 )
  {
    G4cout << "Born excitation model is initialized " << G4endl
    << "Energy range: "
    << LowEnergyLimit() / eV << " eV - "
    << HighEnergyLimit() / keV << " keV for "
    << particle->GetParticleName()
    << G4endl;
  }

  // Initialize water density pointer
  fpMolWaterDensity = G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));

  if (isInitialised)
  { return;}
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNABornExcitationModel::CrossSectionPerVolume(const G4Material* material,
                                                         const G4ParticleDefinition* particleDefinition,
                                                         G4double ekin,
                                                         G4double,
                                                         G4double)
{
  if (verboseLevel > 3)
  G4cout << "Calling CrossSectionPerVolume() of G4DNABornExcitationModel"
  << G4endl;

  if (
      particleDefinition != G4Proton::ProtonDefinition()
      &&
      particleDefinition != G4Electron::ElectronDefinition()
  )

  return 0;

  // Calculate total cross section for model

  G4double lowLim = 0;
  G4double highLim = 0;
  G4double sigma=0;

  G4double waterDensity = (*fpMolWaterDensity)[material->GetIndex()];

  if(waterDensity!= 0.0)
  //  if (material == nistwater || material->GetBaseMaterial() == nistwater)
  {
    const G4String& particleName = particleDefinition->GetParticleName();

    std::map< G4String,G4double,std::less<G4String> >::iterator pos1;
    pos1 = lowEnergyLimit.find(particleName);
    if (pos1 != lowEnergyLimit.end())
    {
      lowLim = pos1->second;
    }

    std::map< G4String,G4double,std::less<G4String> >::iterator pos2;
    pos2 = highEnergyLimit.find(particleName);
    if (pos2 != highEnergyLimit.end())
    {
      highLim = pos2->second;
    }

    if (ekin >= lowLim && ekin < highLim)
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
        G4Exception("G4DNABornExcitationModel::CrossSectionPerVolume","em0002",
            FatalException,"Model not applicable to particle type.");
      }
    }

    if (verboseLevel > 2)
    {
      G4cout << "__________________________________" << G4endl;
      G4cout << "=== G4DNABornExcitationModel - XS INFO START" << G4endl;
      G4cout << "=== Kinetic energy(eV)=" << ekin/eV << " particle : " << particleName << G4endl;
      G4cout << "=== Cross section per water molecule (cm^2)=" << sigma/cm/cm << G4endl;
      G4cout << "=== Cross section per water molecule (cm^-1)=" << sigma*waterDensity/(1./cm) << G4endl;
      //      G4cout << " - Cross section per water molecule (cm^-1)=" << sigma*material->GetAtomicNumDensityVector()[1]/(1./cm) << G4endl;
      G4cout << "=== G4DNABornExcitationModel - XS INFO END" << G4endl;
    }

  } // if (waterMaterial)

  return sigma*waterDensity;
  // return sigma*material->GetAtomicNumDensityVector()[1];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNABornExcitationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
                                                 const G4MaterialCutsCouple* /*couple*/,
                                                 const G4DynamicParticle* aDynamicParticle,
                                                 G4double,
                                                 G4double)
{

  if (verboseLevel > 3)
  G4cout << "Calling SampleSecondaries() of G4DNABornExcitationModel" << G4endl;

  G4double k = aDynamicParticle->GetKineticEnergy();

  const G4String& particleName = aDynamicParticle->GetDefinition()->GetParticleName();

  G4int level = RandomSelect(k,particleName);
  G4double excitationEnergy = waterStructure.ExcitationEnergy(level);
  G4double newEnergy = k - excitationEnergy;

  if (newEnergy > 0)
  {
    fParticleChangeForGamma->ProposeMomentumDirection(aDynamicParticle->GetMomentumDirection());
    fParticleChangeForGamma->SetProposedKineticEnergy(newEnergy);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(excitationEnergy);
  }

  const G4Track * theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
  G4DNAChemistryManager::Instance()->CreateWaterMolecule(eExcitedMolecule,
      level,
      theIncomingTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNABornExcitationModel::RandomSelect(G4double k,
                                             const G4String& particle)
{
  G4int level = 0;

  std::map<G4String, G4DNACrossSectionDataSet*, std::less<G4String> >::iterator pos;
  pos = tableData.find(particle);

  if (pos != tableData.end())
  {
    G4DNACrossSectionDataSet* table = pos->second;

    if (table != 0)
    {
      G4double* valuesBuffer = new G4double[table->NumberOfComponents()];
      const size_t n(table->NumberOfComponents());
      size_t i(n);
      G4double value = 0.;

      while (i > 0)
      {
        i--;
        valuesBuffer[i] = table->GetComponent(i)->FindValue(k);
        value += valuesBuffer[i];
      }

      value *= G4UniformRand();

      i = n;

      while (i > 0)
      {
        i--;

        if (valuesBuffer[i] > value)
        {
          delete[] valuesBuffer;
          return i;
        }
        value -= valuesBuffer[i];
      }

      if (valuesBuffer) delete[] valuesBuffer;

    }
  }
  else
  {
    G4Exception("G4DNABornExcitationModel::RandomSelect", "em0002",
                FatalException, "Model not applicable to particle type.");
  }
  return level;
}

