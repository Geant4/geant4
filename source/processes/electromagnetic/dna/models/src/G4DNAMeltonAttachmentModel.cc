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
// $Id: G4DNAMeltonAttachmentModel.cc 98733 2016-08-09 10:51:58Z gcosmo $
//

// Created by Z. Francis

#include "G4DNAMeltonAttachmentModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//#define MELTON_VERBOSE // prevent checking conditions at run time

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAMeltonAttachmentModel::G4DNAMeltonAttachmentModel(const G4ParticleDefinition*,
                                                       const G4String& nam) :
    G4VEmModel(nam), isInitialised(false)
{
  fpWaterDensity = 0;

  SetLowEnergyLimit(4.*eV);
  SetHighEnergyLimit(13.*eV);

  verboseLevel = 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

#ifdef MELTON_VERBOSE
  if (verboseLevel > 0)
  {
    G4cout << "Melton Attachment model is constructed "
           << G4endl
           << "Energy range: "
           << LowEnergyLimit() / eV << " eV - "
           << HighEnergyLimit() / eV << " eV"
           << G4endl;
  }
#endif
  
  fParticleChangeForGamma = 0;
  fDissociationFlag = true;
  fData = 0;

  // Selection of stationary mode

  statCode = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAMeltonAttachmentModel::~G4DNAMeltonAttachmentModel()
{
  if(fData) delete fData;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAMeltonAttachmentModel::Initialise(const G4ParticleDefinition* particle,
                                            const G4DataVector& /*cuts*/)
{
#ifdef MELTON_VERBOSE
  if (verboseLevel > 3)
    G4cout
      << "Calling G4DNAMeltonAttachmentModel::Initialise()" << G4endl;
#endif

  // ONLY ELECTRON
  
  if(particle->GetParticleName() != "e-")
  {
    G4Exception("G4DNAMeltonAttachmentModel::Initialise",
                "em0002",
                FatalException,
                "Model not applicable to particle type.");
  }
  
  // Energy limits

  if (LowEnergyLimit() < 4.*eV)
  {
    G4ExceptionDescription errMsg;
    errMsg << "G4DNAMeltonAttachmentModel: low energy limit increased from " <<
    LowEnergyLimit()/eV << " eV to " << 4.  << " eV" << G4endl;
    
    G4Exception("G4DNAMeltonAttachmentModel::Initialise",
                "Melton_LowerEBoundary",
                JustWarning,
                errMsg);
    
    SetLowEnergyLimit(4*eV);
  }

  if (HighEnergyLimit() > 13.*eV)
  {
    G4ExceptionDescription errMsg;
    errMsg << "G4DNAMeltonAttachmentModel: high energy limit decreased from " <<
    HighEnergyLimit()/eV << " eV to " << 13. << " eV" << G4endl;
    
    G4Exception("G4DNAMeltonAttachmentModel::Initialise",
                "Melton_HigherEBoundary",
                JustWarning,
                errMsg);
    
    SetHighEnergyLimit(13.*eV);
  }

  // Reading of data files

  G4double scaleFactor = 1e-18*cm2;

  // For total cross section
  G4String fileElectron("dna/sigma_attachment_e_melton");

  fData = new G4DNACrossSectionDataSet(new G4LogLogInterpolation(),
                                     eV, scaleFactor);
  fData->LoadData(fileElectron);


#ifdef MELTON_VERBOSE
  if( verboseLevel >0)
  {
    if (verboseLevel > 2)
    {
      G4cout << "Loaded cross section data for Melton Attachment model" << G4endl;
    }
    
    G4cout << "Melton Attachment model is initialized " << G4endl
    << "Energy range: "
    << LowEnergyLimit() / eV << " eV - "
    << HighEnergyLimit() / eV << " eV"
    << G4endl;
  }
#endif 
  
  // Initialize water density pointer
  fpWaterDensity = G4DNAMolecularMaterial::Instance()->
      GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));

  if (isInitialised)
  {
    return;
  }
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4DNAMeltonAttachmentModel::CrossSectionPerVolume(const G4Material* material,
                                                const G4ParticleDefinition*,
                                                G4double ekin,
                                                G4double,
                                                G4double)
{
#ifdef MELTON_VERBOSE
  if (verboseLevel > 3)
    G4cout
      << "Calling CrossSectionPerVolume() of G4DNAMeltonAttachmentModel"
      << G4endl;
#endif

  // Calculate total cross section for model

  G4double sigma = 0.;

  G4double waterDensity = (*fpWaterDensity)[material->GetIndex()];

  if(waterDensity != 0.0)
  {
    if (ekin >= LowEnergyLimit() && ekin < HighEnergyLimit())
      // necessaire ?
    {
      sigma = fData->FindValue(ekin);
    }

#ifdef MELTON_VERBOSE
    if (verboseLevel > 2)
    {
      G4cout << "__________________________________" << G4endl;
      G4cout << "=== G4DNAMeltonAttachmentModel - XS INFO START" << G4endl;
      G4cout << "--- Kinetic energy(eV)=" << ekin/eV
            << " particle : " << particleDefinition->GetParticleName()
            << G4endl;
      G4cout << "--- Cross section per water molecule (cm^2)="
             << sigma/cm/cm << G4endl;
      G4cout << "--- Cross section per water molecule (cm^-1)="
             << sigma*waterDensity/(1./cm) << G4endl;
      G4cout << "--- G4DNAMeltonAttachmentModel - XS INFO END" << G4endl;
    }
#endif
  } // if water

  return sigma*waterDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void
G4DNAMeltonAttachmentModel::
SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
                  const G4MaterialCutsCouple* /*couple*/,
                  const G4DynamicParticle* aDynamicElectron,
                  G4double,
                  G4double)
{
  
#ifdef MELTON_VERBOSE
  if (verboseLevel > 3)
    G4cout
    << "Calling SampleSecondaries() of G4DNAMeltonAttachmentModel" << G4endl;
#endif
  
  // Electron is killed
  
  G4double electronEnergy0 = aDynamicElectron->GetKineticEnergy();

  if (!statCode)     
  {     
      fParticleChangeForGamma->SetProposedKineticEnergy(0.);
      fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(electronEnergy0);
  }

  else 
  {
      fParticleChangeForGamma->SetProposedKineticEnergy(electronEnergy0);
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(electronEnergy0);
  }
  
  if(fDissociationFlag)
  {
    G4DNAChemistryManager::Instance()->
      CreateWaterMolecule(eDissociativeAttachment,
                          -1,
                          fParticleChangeForGamma->GetCurrentTrack());
  }
  return;
}
