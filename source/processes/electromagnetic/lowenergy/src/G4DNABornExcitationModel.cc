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
// $Id: G4DNABornExcitationModel.cc,v 1.7 2009/08/31 14:03:29 sincerti Exp $
// GEANT4 tag $Name: geant4-09-03 $
//

#include "G4DNABornExcitationModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNABornExcitationModel::G4DNABornExcitationModel(const G4ParticleDefinition*,
                                             const G4String& nam)
:G4VEmModel(nam),isInitialised(false)
{

  lowEnergyLimit = 500 * keV; 
  highEnergyLimit = 100 * MeV;
  SetLowEnergyLimit(lowEnergyLimit);
  SetHighEnergyLimit(highEnergyLimit);

  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods
  
  //
  
  table = 0;
  
  //
  
  if( verboseLevel>0 ) 
  { 
    G4cout << "Born excitation model is constructed " << G4endl
           << "Energy range: "
           << lowEnergyLimit / keV << " keV - "
           << highEnergyLimit / MeV << " MeV"
           << G4endl;
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNABornExcitationModel::~G4DNABornExcitationModel()
{ 
  // Cross section
  delete table;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNABornExcitationModel::Initialise(const G4ParticleDefinition* /*particle*/,
                                       const G4DataVector& /*cuts*/)
{

  if (verboseLevel > 3)
    G4cout << "Calling G4DNABornExcitationModel::Initialise()" << G4endl;

  // Energy limits
  
  if (LowEnergyLimit() < lowEnergyLimit)
  {
    G4cout << "G4DNABornExcitationModel: low energy limit increased from " << 
	LowEnergyLimit()/keV << " keV to " << lowEnergyLimit/keV << " keV" << G4endl;
    SetLowEnergyLimit(lowEnergyLimit);
  }

  if (HighEnergyLimit() > highEnergyLimit)
  {
    G4cout << "G4DNABornExcitationModel: high energy limit decreased from " << 
        HighEnergyLimit()/MeV << " MeV to " << highEnergyLimit/MeV << " MeV" << G4endl;
    SetHighEnergyLimit(highEnergyLimit);
  }

  //
  
  if (table == 0)
  {
    table = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,(1e-22/3.343)*m*m );
    table->LoadData("dna/sigma_excitation_p_born");
  }

  if( verboseLevel>0 ) 
  { 
    G4cout << "Born excitation model is initialized " << G4endl
           << "Energy range: "
           << LowEnergyLimit() / keV << " keV - "
           << HighEnergyLimit() / MeV << " MeV " << G4endl;
  }
  
  if(!isInitialised) 
  {
    isInitialised = true;
  
    if(pParticleChange)
      fParticleChangeForGamma = reinterpret_cast<G4ParticleChangeForGamma*>(pParticleChange);
    else
      fParticleChangeForGamma = new G4ParticleChangeForGamma();
  }    

  // InitialiseElementSelectors(particle,cuts);
  
  // Test if water material

  flagMaterialIsWater= false;
  densityWater = 0;

  const G4ProductionCutsTable* theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();

  if(theCoupleTable) 
  {
    G4int numOfCouples = theCoupleTable->GetTableSize();
  
    if(numOfCouples>0) 
    {
	  for (G4int i=0; i<numOfCouples; i++) 
	  {
	    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
	    const G4Material* material = couple->GetMaterial();

            if (material->GetName() == "G4_WATER") 
            {
              G4double density = material->GetAtomicNumDensityVector()[1];
	      flagMaterialIsWater = true; 
	      densityWater = density; 
	      
	      if (verboseLevel > 3) 
              G4cout << "****** Water material is found with density(cm^-3)=" << density/(cm*cm*cm) << G4endl;
            }
  
          }

    } // if(numOfCouples>0)

  } // if (theCoupleTable)

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNABornExcitationModel::CrossSectionPerVolume(const G4Material*,
					   const G4ParticleDefinition* particleDefinition,
					   G4double k,
					   G4double,
					   G4double)
{
  if (verboseLevel > 3)
    G4cout << "Calling CrossSectionPerVolume() of G4DNABornExcitationModel" << G4endl;

 // Calculate total cross section for model

  G4double crossSection=0;
  
  if (flagMaterialIsWater)
  {
    if (particleDefinition == G4Proton::ProtonDefinition())
    {
      if (k >= lowEnergyLimit && k < highEnergyLimit)
      {
        crossSection = table->FindValue(k); 
      }

      if (verboseLevel > 3)
      {
        G4cout << "---> Kinetic energy(keV)=" << k/keV << G4endl;
        G4cout << " - Cross section per water molecule (cm^2)=" << crossSection/cm/cm << G4endl;
        G4cout << " - Cross section per water molecule (cm^-1)=" << crossSection*densityWater/(1./cm) << G4endl;
      } 
    }
  } // if (flagMaterialIsWater)

 return crossSection*densityWater;		   

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
  
  G4int level = RandomSelect(k);
  G4double excitationEnergy = waterStructure.ExcitationEnergy(level);
  G4double newEnergy = k - excitationEnergy;
  
  if (newEnergy > 0)
  {
      fParticleChangeForGamma->ProposeMomentumDirection(aDynamicParticle->GetMomentumDirection());
      fParticleChangeForGamma->SetProposedKineticEnergy(newEnergy);
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(excitationEnergy);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNABornExcitationModel::RandomSelect(G4double k)
{   
  G4int level = 0;

  G4double* valuesBuffer = new G4double[table->NumberOfComponents()];

  const size_t n(table->NumberOfComponents());
  size_t i(n);
  G4double value = 0.;
  
  while (i>0)
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

  return level;
}


