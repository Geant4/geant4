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
// $Id: G4DNAMeltonAttachmentModel.cc,v 1.2 2010-09-15 05:47:33 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// Created by Z. Francis

#include "G4DNAMeltonAttachmentModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAMeltonAttachmentModel::G4DNAMeltonAttachmentModel(const G4ParticleDefinition*,
                                             const G4String& nam)
:G4VEmModel(nam),isInitialised(false)
{

  lowEnergyLimit = 4 * eV; 
  lowEnergyLimitOfModel = 4 * eV; 
  highEnergyLimit = 13 * eV;
  SetLowEnergyLimit(lowEnergyLimit);
  SetHighEnergyLimit(highEnergyLimit);

  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods
  
  if( verboseLevel>0 ) 
  { 
    G4cout << "Melton Attachment model is constructed " << G4endl
           << "Energy range: "
           << lowEnergyLimit / eV << " eV - "
           << highEnergyLimit / eV << " eV"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAMeltonAttachmentModel::~G4DNAMeltonAttachmentModel()
{  
  // For total cross section
  
  std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
  
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

  if (verboseLevel > 3)
    G4cout << "Calling G4DNAMeltonAttachmentModel::Initialise()" << G4endl;

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
 
  if (electronDef != 0)
  {
    // For total cross section
    
    electron = electronDef->GetParticleName();

    tableFile[electron] = fileElectron;

    G4DNACrossSectionDataSet* tableE = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
    tableE->LoadData(fileElectron);
    tableData[electron] = tableE;
    
  }
  else G4Exception("G4DNAMeltonAttachmentModel::Initialise: electron is not defined");
  
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

  if(!isInitialised) 
  {
    isInitialised = true;
  
    if(pParticleChange)
      fParticleChangeForGamma = reinterpret_cast<G4ParticleChangeForGamma*>(pParticleChange);
    else
      fParticleChangeForGamma = new G4ParticleChangeForGamma();
  }    

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAMeltonAttachmentModel::CrossSectionPerVolume(const G4Material* material,
					   const G4ParticleDefinition* p,
					   G4double ekin,
					   G4double,
					   G4double)
{
  if (verboseLevel > 3)
    G4cout << "Calling CrossSectionPerVolume() of G4DNAMeltonAttachmentModel" << G4endl;

 // Calculate total cross section for model

 G4double sigma=0;
 
 if (material->GetName() == "G4_WATER")
 {
  const G4String& particleName = p->GetParticleName();

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
	    G4Exception("G4DNAMeltonAttachmentModel::ComputeCrossSectionPerVolume: attempting to calculate cross section for wrong particle");
	}
  }

  if (verboseLevel > 3)
  {
    G4cout << "---> Kinetic energy(eV)=" << ekin/eV << G4endl;
    G4cout << " - Cross section per water molecule (cm^2)=" << sigma/cm/cm << G4endl;
    G4cout << " - Cross section per water molecule (cm^-1)=" << sigma*material->GetAtomicNumDensityVector()[1]/(1./cm) << G4endl;
  } 

 } // if water
        
 return sigma*material->GetAtomicNumDensityVector()[1];		   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAMeltonAttachmentModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
					      const G4MaterialCutsCouple* /*couple*/,
					      const G4DynamicParticle* aDynamicElectron,
					      G4double,
					      G4double)
{

  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4DNAMeltonAttachmentModel" << G4endl;

  // Electron is killed
  
  G4double electronEnergy0 = aDynamicElectron->GetKineticEnergy();
  fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
  fParticleChangeForGamma->ProposeLocalEnergyDeposit(electronEnergy0);
  
  return ;
}








