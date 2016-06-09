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
// $Id: G4LivermoreRayleighModel.cc,v 1.9 2010-12-27 17:45:12 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Sebastien Inserti
//         30 October 2008
//
// History:
// --------
// 18 Apr 2009   V Ivanchenko Cleanup initialisation and generation of secondaries:
//                  - apply internal high-energy limit only in constructor 
//                  - do not apply low-energy limit (default is 0)
//                  - remove GetMeanFreePath method and table
//                  - remove initialisation of element selector 
//                  - use G4ElementSelector
// 26 Dec 2010   V Ivanchenko Load data tables only once to avoid memory leak

#include "G4LivermoreRayleighModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermoreRayleighModel::G4LivermoreRayleighModel(const G4ParticleDefinition*,
						   const G4String& nam)
  :G4VEmModel(nam),isInitialised(false),meanFreePathTable(0),
   formFactorData(0),crossSectionHandler(0)
{
  lowEnergyLimit = 250 * eV; 
  highEnergyLimit = 100 * GeV;
  
  //  SetLowEnergyLimit(lowEnergyLimit);
  SetHighEnergyLimit(highEnergyLimit);
  //
  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if(verboseLevel > 0) {
    G4cout << "Livermore Rayleigh is constructed " << G4endl
	   << "Energy range: "
	   << lowEnergyLimit / eV << " eV - "
	   << highEnergyLimit / GeV << " GeV"
	   << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermoreRayleighModel::~G4LivermoreRayleighModel()
{  
  delete crossSectionHandler;
  delete formFactorData;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreRayleighModel::Initialise(const G4ParticleDefinition* particle,
					  const G4DataVector& cuts)
{
  if (verboseLevel > 3) {
    G4cout << "Calling G4LivermoreRayleighModel::Initialise()" << G4endl;
  }

  if (crossSectionHandler)
  {
    crossSectionHandler->Clear();
    delete crossSectionHandler;
  }
  delete formFactorData;

  // Data are read for all materials
  
  crossSectionHandler = new G4CrossSectionHandler;
  //  crossSectionHandler->Clear();
  G4String crossSectionFile = "rayl/re-cs-";
  crossSectionHandler->LoadData(crossSectionFile);

  G4VDataSetAlgorithm* ffInterpolation = new G4LogLogInterpolation;
  G4String formFactorFile = "rayl/re-ff-";
  formFactorData = new G4CompositeEMDataSet(ffInterpolation,1.,1.);
  formFactorData->LoadData(formFactorFile);
  
  InitialiseElementSelectors(particle,cuts);

  //  
  if (verboseLevel > 2) {
    G4cout << "Loaded cross section files for Livermore Rayleigh model" << G4endl;
  }
  if (verboseLevel > 0) { 
    G4cout << "Livermore Rayleigh model is initialized " << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / eV << " eV - "
	   << HighEnergyLimit() / GeV << " GeV"
	   << G4endl;
  }

  if(isInitialised) { return; }
  fParticleChange = GetParticleChangeForGamma();
  isInitialised = true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermoreRayleighModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double GammaEnergy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
  if (verboseLevel > 3) {
    G4cout << "Calling CrossSectionPerAtom() of G4LivermoreRayleighModel" << G4endl;
  }
  if (GammaEnergy < lowEnergyLimit || GammaEnergy > highEnergyLimit) {
    return 0.0;
  }
  G4double cs = crossSectionHandler->FindValue(G4int(Z), GammaEnergy);
  return cs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreRayleighModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
					      const G4MaterialCutsCouple* couple,
					      const G4DynamicParticle* aDynamicGamma,
					      G4double,
					      G4double)
{
  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4LivermoreRayleighModel" << G4endl;

  G4double photonEnergy0 = aDynamicGamma->GetKineticEnergy();

  // absorption of low-energy gamma  
  if (photonEnergy0 <= lowEnergyLimit)
    {
      fParticleChange->ProposeTrackStatus(fStopAndKill);
      fParticleChange->SetProposedKineticEnergy(0.);
      fParticleChange->ProposeLocalEnergyDeposit(photonEnergy0);
      return ;
    }

  G4ParticleMomentum photonDirection0 = aDynamicGamma->GetMomentumDirection();

  // Select randomly one element in the current material
  //  G4int Z = crossSectionHandler->SelectRandomAtom(couple,photonEnergy0);
  const G4ParticleDefinition* particle =  aDynamicGamma->GetDefinition();
  const G4Element* elm = SelectRandomAtom(couple,particle,photonEnergy0);
  G4int Z = (G4int)elm->GetZ();

  // Sample the angle of the scattered photon

  G4double wlPhoton = h_Planck*c_light/photonEnergy0;

  G4double gReject,x,dataFormFactor;
  G4double randomFormFactor;
  G4double cosTheta;
  G4double sinTheta;
  G4double fcostheta;

  do
    {
      do
	{
	  cosTheta = 2. * G4UniformRand() - 1.;
	  fcostheta = ( 1. + cosTheta*cosTheta)/2.;
	} while (fcostheta < G4UniformRand());

      if (photonEnergy0 > 5)
	{
	  cosTheta = 1.;
	}

      G4double sinThetaHalf = std::sqrt((1. - cosTheta) / 2.);
      x = sinThetaHalf / (wlPhoton/cm);
      if (x > 1.e+005)
	{
	  dataFormFactor = formFactorData->FindValue(x,Z-1);
	}
      else
	{
	  dataFormFactor = formFactorData->FindValue(0.,Z-1);
	}
      randomFormFactor = G4UniformRand() * Z * Z;
      sinTheta = std::sqrt(1. - cosTheta*cosTheta);
      gReject = dataFormFactor * dataFormFactor;

    } while( gReject < randomFormFactor);

  // Scattered photon angles. ( Z - axis along the parent photon)
  G4double phi = twopi * G4UniformRand() ;
  G4double dirX = sinTheta*std::cos(phi);
  G4double dirY = sinTheta*std::sin(phi);
  G4double dirZ = cosTheta;

  // Update G4VParticleChange for the scattered photon
  G4ThreeVector photonDirection1(dirX, dirY, dirZ);
  photonDirection1.rotateUz(photonDirection0);
  fParticleChange->ProposeMomentumDirection(photonDirection1);

  fParticleChange->SetProposedKineticEnergy(photonEnergy0); 
}


