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
// Author: Luciano Pandola
//         on base of G4LowEnergyBremsstrahlung developed by A.Forti and V.Ivanchenko
//
// History:
// --------
// 03 Mar 2009   L Pandola    Migration from process to model 
// 12 Apr 2009   V Ivanchenko Cleanup initialisation and generation of secondaries:
//                  - apply internal high-energy limit only in constructor 
//                  - do not apply low-energy limit (default is 0)
//                  - added MinEnergyCut method
//                  - do not change track status
//                  - do not initialize element selectors
//                  - use cut value from the interface 
//                  - fixed bug in sampling of angles between keV and MeV
// 19 May 2009   L Pandola    Explicitely set to zero pointers deleted in 
//                            Initialise(), since they might be checked later on
//

#include "G4LivermoreBremsstrahlungModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"

#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4SemiLogInterpolation.hh"
//
#include "G4VEmAngularDistribution.hh"
#include "G4ModifiedTsai.hh"
#include "G4Generator2BS.hh"
//#include "G4Generator2BN.hh"
//
#include "G4BremsstrahlungCrossSectionHandler.hh"
//
#include "G4VEnergySpectrum.hh"
#include "G4eBremsstrahlungSpectrum.hh"
#include "G4VEMDataSet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4LivermoreBremsstrahlungModel::G4LivermoreBremsstrahlungModel(const G4ParticleDefinition*,
							       const G4String& nam)
  :G4VEmModel(nam),fParticleChange(0),isInitialised(false),
   crossSectionHandler(0),energySpectrum(0)
{
  fIntrinsicLowEnergyLimit = 10.0*eV;
  fIntrinsicHighEnergyLimit = 100.0*GeV;
  fNBinEnergyLoss = 360;
  //  SetLowEnergyLimit(fIntrinsicLowEnergyLimit);
  SetHighEnergyLimit(fIntrinsicHighEnergyLimit);
  //
  verboseLevel = 0;
  SetAngularDistribution(new G4Generator2BS());
  //
  //generatorName = "tsai";
  //angularDistribution = new G4ModifiedTsai("TsaiGenerator"); //default generator
  //
  //TsaiAngularDistribution = new G4ModifiedTsai("TsaiGenerator");
  //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermoreBremsstrahlungModel::~G4LivermoreBremsstrahlungModel()
{
  if (crossSectionHandler) delete crossSectionHandler;
  if (energySpectrum) delete energySpectrum;
  energyBins.clear();
  //delete angularDistribution;
  //delete TsaiAngularDistribution;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreBremsstrahlungModel::Initialise(const G4ParticleDefinition* particle,
						const G4DataVector& cuts)
{
  //Check that the Livermore Bremsstrahlung is NOT attached to e+
  if (particle != G4Electron::Electron())
    {
      G4Exception("G4LivermoreBremsstrahlungModel::Initialise",
		    "em0002",FatalException,"Livermore Bremsstrahlung Model is applicable only to electrons");
    }
  //Prepare energy spectrum
  if (energySpectrum) 
    {
      delete energySpectrum;
      energySpectrum = 0;
    }

  energyBins.clear();
  for(size_t i=0; i<15; i++) 
    {
      G4double x = 0.1*((G4double)i);
      if(i == 0)  x = 0.01;
      if(i == 10) x = 0.95;
      if(i == 11) x = 0.97;
      if(i == 12) x = 0.99;
      if(i == 13) x = 0.995;
      if(i == 14) x = 1.0;
      energyBins.push_back(x);
    }
  const G4String dataName("/brem/br-sp.dat");
  energySpectrum = new G4eBremsstrahlungSpectrum(energyBins,dataName);
  
  if (verboseLevel > 0)
    G4cout << "G4eBremsstrahlungSpectrum is initialized" << G4endl;

  //Initialize cross section handler
  if (crossSectionHandler) 
    {
      delete crossSectionHandler;
      crossSectionHandler = 0;
    }
  G4VDataSetAlgorithm* interpolation = 0;//new G4SemiLogInterpolation();
  crossSectionHandler = new G4BremsstrahlungCrossSectionHandler(energySpectrum,interpolation);
  crossSectionHandler->Initialise(0,LowEnergyLimit(),HighEnergyLimit(),
				  fNBinEnergyLoss);
  crossSectionHandler->Clear();
  crossSectionHandler->LoadShellData("brem/br-cs-");
  //This is used to retrieve cross section values later on
  G4VEMDataSet* p = crossSectionHandler->BuildMeanFreePathForMaterials(&cuts);
  delete p;  
 
  if (verboseLevel > 0)
    {
      G4cout << "Livermore Bremsstrahlung model is initialized " << G4endl
	     << "Energy range: "
	     << LowEnergyLimit() / keV << " keV - "
	     << HighEnergyLimit() / GeV << " GeV"
	     << G4endl;
    }

  if (verboseLevel > 1)
    {
      G4cout << "Cross section data: " << G4endl; 
      crossSectionHandler->PrintData();
      G4cout << "Parameters: " << G4endl;
      energySpectrum->PrintData();
    }

  if(isInitialised) return;
  fParticleChange = GetParticleChangeForLoss();
  isInitialised = true; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermoreBremsstrahlungModel::MinEnergyCut(const G4ParticleDefinition*,
						      const G4MaterialCutsCouple*)
{
  return 250.*eV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4LivermoreBremsstrahlungModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
							   G4double energy,
							   G4double Z, G4double,
							   G4double cutEnergy, 
							   G4double)
{
  G4int iZ = (G4int) Z;
  if (!crossSectionHandler)
    {
      G4Exception("G4LivermoreBremsstrahlungModel::ComputeCrossSectionPerAtom",
		    "em1007",FatalException,"The cross section handler is not correctly initialized");
      return 0;
    }
  
  //The cut is already included in the crossSectionHandler
  G4double cs = 
    crossSectionHandler->GetCrossSectionAboveThresholdForElement(energy,cutEnergy,iZ);

  if (verboseLevel > 1)
    {
      G4cout << "G4LivermoreBremsstrahlungModel " << G4endl;
      G4cout << "Cross section for gamma emission > " << cutEnergy/keV << " keV at " <<
	energy/keV << " keV and Z = " << iZ << " --> " << cs/barn << " barn" << G4endl;
    }
  return cs;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermoreBremsstrahlungModel::ComputeDEDXPerVolume(const G4Material* material,
                             		   const G4ParticleDefinition* ,
                               		   G4double kineticEnergy,
                               		   G4double cutEnergy)
{
  G4double sPower = 0.0;

  const G4ElementVector* theElementVector = material->GetElementVector();
  size_t NumberOfElements = material->GetNumberOfElements() ;
  const G4double* theAtomicNumDensityVector =
                    material->GetAtomicNumDensityVector();

  // loop for elements in the material
  for (size_t iel=0; iel<NumberOfElements; iel++ ) 
    {
      G4int iZ = (G4int)((*theElementVector)[iel]->GetZ());
      G4double e = energySpectrum->AverageEnergy(iZ, 0.0,cutEnergy,
						 kineticEnergy);
      G4double cs= crossSectionHandler->FindValue(iZ,kineticEnergy);
      sPower   += e * cs * theAtomicNumDensityVector[iel];
    }

  if (verboseLevel > 2)
    {
      G4cout << "G4LivermoreBremsstrahlungModel " << G4endl;
      G4cout << "Stopping power < " << cutEnergy/keV << " keV at " << 
	kineticEnergy/keV << " keV = " << sPower/(keV/mm) << " keV/mm" << G4endl;
    }
    
  return sPower;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreBremsstrahlungModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
					      const G4MaterialCutsCouple* couple,
					      const G4DynamicParticle* aDynamicParticle,
					      G4double energyCut,
					      G4double)
{
  
  G4double kineticEnergy = aDynamicParticle->GetKineticEnergy();

  // this is neede for pathalogical cases of no ionisation
  if (kineticEnergy <= fIntrinsicLowEnergyLimit)
    {
      fParticleChange->SetProposedKineticEnergy(0.);
      fParticleChange->ProposeLocalEnergyDeposit(kineticEnergy);
      return;
    }

  //Sample material
  G4int Z = crossSectionHandler->SelectRandomAtom(couple, kineticEnergy);

  //Sample gamma energy
  G4double tGamma = energySpectrum->SampleEnergy(Z, energyCut, kineticEnergy, kineticEnergy);
  //nothing happens
  if (tGamma == 0.) { return; }

  G4double totalEnergy = kineticEnergy + electron_mass_c2;
  G4double finalEnergy = kineticEnergy - tGamma; // electron final energy  

  //Sample gamma direction
  G4ThreeVector gammaDirection = 
    GetAngularDistribution()->SampleDirection(aDynamicParticle, 
					      totalEnergy-tGamma,
					      Z, 
					      couple->GetMaterial());

  G4ThreeVector electronDirection = aDynamicParticle->GetMomentumDirection();

  //Update the incident particle    
  if (finalEnergy < 0.) 
    {
      // Kinematic problem
      tGamma = kineticEnergy;
      fParticleChange->SetProposedKineticEnergy(0.);
    }
  else
    {
      G4double momentum = std::sqrt((totalEnergy + electron_mass_c2)*kineticEnergy);
      G4double finalX = momentum*electronDirection.x() - tGamma*gammaDirection.x();
      G4double finalY = momentum*electronDirection.y() - tGamma*gammaDirection.y();
      G4double finalZ = momentum*electronDirection.z() - tGamma*gammaDirection.z();
      G4double norm = 1./std::sqrt(finalX*finalX + finalY*finalY + finalZ*finalZ);
      
      fParticleChange->ProposeMomentumDirection(finalX*norm, finalY*norm, finalZ*norm);
      fParticleChange->SetProposedKineticEnergy(finalEnergy);
    }

  //Generate the bremsstrahlung gamma
  G4DynamicParticle* aGamma= new G4DynamicParticle (G4Gamma::Gamma(),
						    gammaDirection, tGamma);
  fvect->push_back(aGamma);

  if (verboseLevel > 1)
    {
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Energy balance from G4LivermoreBremsstrahlung" << G4endl;
      G4cout << "Incoming primary energy: " << kineticEnergy/keV << " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Outgoing primary energy: " << finalEnergy/keV << " keV" << G4endl;
      G4cout << "Gamma ray " << tGamma/keV << " keV" << G4endl;
      G4cout << "Total final state: " << (finalEnergy+tGamma)/keV << " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
    }
  if (verboseLevel > 0)
    {
      G4double energyDiff = std::fabs(finalEnergy+tGamma-kineticEnergy);
      if (energyDiff > 0.05*keV)
	G4cout << "G4LivermoreBremsstrahlung WARNING: problem with energy conservation: " 
	       << (finalEnergy+tGamma)/keV << " keV (final) vs. " 
	       << kineticEnergy/keV << " keV (initial)" << G4endl;
    }
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
/*
void 
G4LivermoreBremsstrahlungModel::SetAngularGenerator(G4VBremAngularDistribution* distribution)
{
  if(angularDistribution == distribution) return;
  if(angularDistribution) delete angularDistribution;
  angularDistribution = distribution;
  angularDistribution->PrintGeneratorInformation();
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 /*
void G4LivermoreBremsstrahlungModel::SetAngularGenerator(const G4String& theGenName)
{
  if(theGenName == generatorName) return;
  if (theGenName == "tsai") 
    {
      delete angularDistribution;
      angularDistribution = new G4ModifiedTsai("TsaiGenerator");
      generatorName = theGenName;
    }
  else if (theGenName == "2bn")
    {
      delete angularDistribution;
      angularDistribution = new G4Generator2BN("2BNGenerator");
      generatorName = theGenName;
    }
  else if (theGenName == "2bs")
    {
      delete angularDistribution;
      angularDistribution = new G4Generator2BS("2BSGenerator");
      generatorName = theGenName;
    }
  else
    {
      G4cout << "### G4LivermoreBremsstrahlungModel::SetAngularGenerator WARNING:"
	     << " generator <" << theGenName << "> is not known" << G4endl;
      return; 

    }

  angularDistribution->PrintGeneratorInformation();
}
 */
