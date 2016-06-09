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
// $Id: G4LivermoreComptonModel.cc,v 1.8 2010-12-27 17:45:12 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
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
//                  - added protection against numerical problem in energy sampling 
//                  - use G4ElementSelector
// 26 Dec 2010   V Ivanchenko Load data tables only once to avoid memory leak

#include "G4LivermoreComptonModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermoreComptonModel::G4LivermoreComptonModel(const G4ParticleDefinition*,
						 const G4String& nam)
  :G4VEmModel(nam),isInitialised(false),meanFreePathTable(0),
   scatterFunctionData(0),crossSectionHandler(0)
{
  lowEnergyLimit = 250 * eV; 
  highEnergyLimit = 100 * GeV;
  //  SetLowEnergyLimit(lowEnergyLimit);
  SetHighEnergyLimit(highEnergyLimit);

  verboseLevel=0 ;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if(  verboseLevel>0 ) { 
    G4cout << "Livermore Compton model is constructed " << G4endl
	   << "Energy range: "
	   << lowEnergyLimit / eV << " eV - "
	   << highEnergyLimit / GeV << " GeV"
	   << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermoreComptonModel::~G4LivermoreComptonModel()
{  
  delete crossSectionHandler;
  delete scatterFunctionData;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreComptonModel::Initialise(const G4ParticleDefinition* particle,
					 const G4DataVector& cuts)
{
  if (verboseLevel > 3) {
    G4cout << "Calling G4LivermoreComptonModel::Initialise()" << G4endl;
  }

  if (crossSectionHandler)
  {
    crossSectionHandler->Clear();
    delete crossSectionHandler;
  }
  delete scatterFunctionData;

  // Reading of data files - all materials are read
  
  crossSectionHandler = new G4CrossSectionHandler;
  //  crossSectionHandler->Clear();
  G4String crossSectionFile = "comp/ce-cs-";
  crossSectionHandler->LoadData(crossSectionFile);

  G4VDataSetAlgorithm* scatterInterpolation = new G4LogLogInterpolation;
  G4String scatterFile = "comp/ce-sf-";
  scatterFunctionData = new G4CompositeEMDataSet(scatterInterpolation, 1., 1.);
  scatterFunctionData->LoadData(scatterFile);

  // For Doppler broadening
  if(!isInitialised) {
    shellData.SetOccupancyData();
    G4String file = "/doppler/shell-doppler";
    shellData.LoadData(file);
  }

  if (verboseLevel > 2) {
    G4cout << "Loaded cross section files for Livermore Compton model" << G4endl;
  }
 
  InitialiseElementSelectors(particle,cuts);

  if(  verboseLevel>0 ) { 
    G4cout << "Livermore Compton model is initialized " << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / eV << " eV - "
	   << HighEnergyLimit() / GeV << " GeV"
	   << G4endl;
  }
  //  
  if(isInitialised) { return; }
  fParticleChange = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermoreComptonModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double GammaEnergy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
  if (verboseLevel > 3) {
    G4cout << "Calling ComputeCrossSectionPerAtom() of G4LivermoreComptonModel" << G4endl;
  }
  if (GammaEnergy < lowEnergyLimit || GammaEnergy > highEnergyLimit) { return 0.0; }
    
  G4double cs = crossSectionHandler->FindValue(G4int(Z), GammaEnergy);  
  return cs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreComptonModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
						const G4MaterialCutsCouple* couple,
						const G4DynamicParticle* aDynamicGamma,
						G4double, G4double)
{

  // The scattered gamma energy is sampled according to Klein - Nishina formula.
  // then accepted or rejected depending on the Scattering Function multiplied
  // by factor from Klein - Nishina formula.
  // Expression of the angular distribution as Klein Nishina
  // angular and energy distribution and Scattering fuctions is taken from
  // D. E. Cullen "A simple model of photon transport" Nucl. Instr. Meth.
  // Phys. Res. B 101 (1995). Method of sampling with form factors is different
  // data are interpolated while in the article they are fitted.
  // Reference to the article is from J. Stepanek New Photon, Positron
  // and Electron Interaction Data for GEANT in Energy Range from 1 eV to 10
  // TeV (draft).
  // The random number techniques of Butcher & Messel are used
  // (Nucl Phys 20(1960),15).

  G4double photonEnergy0 = aDynamicGamma->GetKineticEnergy();

  if (verboseLevel > 3) {
    G4cout << "G4LivermoreComptonModel::SampleSecondaries() E(MeV)= " 
	   << photonEnergy0/MeV << " in " << couple->GetMaterial()->GetName() 
	   << G4endl;
  }
  
  // low-energy gamma is absorpted by this process
  if (photonEnergy0 <= lowEnergyLimit) 
    {
      fParticleChange->ProposeTrackStatus(fStopAndKill);
      fParticleChange->SetProposedKineticEnergy(0.);
      fParticleChange->ProposeLocalEnergyDeposit(photonEnergy0);
      return ;
    }

  G4double e0m = photonEnergy0 / electron_mass_c2 ;
  G4ParticleMomentum photonDirection0 = aDynamicGamma->GetMomentumDirection();

  // Select randomly one element in the current material
  //  G4int Z = crossSectionHandler->SelectRandomAtom(couple,photonEnergy0);
  const G4ParticleDefinition* particle =  aDynamicGamma->GetDefinition();
  const G4Element* elm = SelectRandomAtom(couple,particle,photonEnergy0);
  G4int Z = (G4int)elm->GetZ();

  G4double epsilon0 = 1. / (1. + 2. * e0m);
  G4double epsilon0Sq = epsilon0 * epsilon0;
  G4double alpha1 = -std::log(epsilon0);
  G4double alpha2 = 0.5 * (1. - epsilon0Sq);

  G4double wlPhoton = h_Planck*c_light/photonEnergy0;

  // Sample the energy of the scattered photon
  G4double epsilon;
  G4double epsilonSq;
  G4double oneCosT;
  G4double sinT2;
  G4double gReject;
  
  do
  {
      if ( alpha1/(alpha1+alpha2) > G4UniformRand())
      {
	// std::pow(epsilon0,G4UniformRand())
        epsilon = std::exp(-alpha1 * G4UniformRand());  
        epsilonSq = epsilon * epsilon;
      }
      else
      {
        epsilonSq = epsilon0Sq + (1. - epsilon0Sq) * G4UniformRand();
        epsilon = std::sqrt(epsilonSq);
      }

      oneCosT = (1. - epsilon) / ( epsilon * e0m);
      sinT2 = oneCosT * (2. - oneCosT);
      G4double x = std::sqrt(oneCosT/2.) / (wlPhoton/cm);
      G4double scatteringFunction = scatterFunctionData->FindValue(x,Z-1);
      gReject = (1. - epsilon * sinT2 / (1. + epsilonSq)) * scatteringFunction;

  } while(gReject < G4UniformRand()*Z);

  G4double cosTheta = 1. - oneCosT;
  G4double sinTheta = std::sqrt (sinT2);
  G4double phi = twopi * G4UniformRand() ;
  G4double dirx = sinTheta * std::cos(phi);
  G4double diry = sinTheta * std::sin(phi);
  G4double dirz = cosTheta ;

  // Doppler broadening -  Method based on:
  // Y. Namito, S. Ban and H. Hirayama, 
  // "Implementation of the Doppler Broadening of a Compton-Scattered Photon 
  // into the EGS4 Code", NIM A 349, pp. 489-494, 1994
  
  // Maximum number of sampling iterations
  G4int maxDopplerIterations = 1000;
  G4double bindingE = 0.;
  G4double photonEoriginal = epsilon * photonEnergy0;
  G4double photonE = -1.;
  G4int iteration = 0;
  G4double eMax = photonEnergy0;
  do
    {
      iteration++;
      // Select shell based on shell occupancy
      G4int shell = shellData.SelectRandomShell(Z);
      bindingE = shellData.BindingEnergy(Z,shell);
      
      eMax = photonEnergy0 - bindingE;
     
      // Randomly sample bound electron momentum 
      // (memento: the data set is in Atomic Units)
      G4double pSample = profileData.RandomSelectMomentum(Z,shell);
      // Rescale from atomic units
      G4double pDoppler = pSample * fine_structure_const;
      G4double pDoppler2 = pDoppler * pDoppler;
      G4double var2 = 1. + oneCosT * e0m;
      G4double var3 = var2*var2 - pDoppler2;
      G4double var4 = var2 - pDoppler2 * cosTheta;
      G4double var = var4*var4 - var3 + pDoppler2 * var3;
      if (var > 0.)
	{
	  G4double varSqrt = std::sqrt(var);        
	  G4double scale = photonEnergy0 / var3;  
          // Random select either root
 	  if (G4UniformRand() < 0.5) photonE = (var4 - varSqrt) * scale;               
	  else photonE = (var4 + varSqrt) * scale;
	} 
      else
	{
	  photonE = -1.;
	}
   } while ( iteration <= maxDopplerIterations && 
	     (photonE < 0. || photonE > eMax || photonE < eMax*G4UniformRand()) );
 
  // End of recalculation of photon energy with Doppler broadening
  // Revert to original if maximum number of iterations threshold has been reached

  if (iteration >= maxDopplerIterations)
    {
      photonE = photonEoriginal;
      bindingE = 0.;
    }

  // Update G4VParticleChange for the scattered photon

  G4ThreeVector photonDirection1(dirx,diry,dirz);
  photonDirection1.rotateUz(photonDirection0);
  fParticleChange->ProposeMomentumDirection(photonDirection1) ;

  G4double photonEnergy1 = photonE;

  if (photonEnergy1 > 0.)
    {
      fParticleChange->SetProposedKineticEnergy(photonEnergy1) ;
    }
  else
    {
      photonEnergy1 = 0.;
      fParticleChange->SetProposedKineticEnergy(0.) ;
      fParticleChange->ProposeTrackStatus(fStopAndKill);   
    }

  // Kinematics of the scattered electron
  G4double eKineticEnergy = photonEnergy0 - photonEnergy1 - bindingE;

  // protection against negative final energy: no e- is created
  if(eKineticEnergy < 0.0) {
    fParticleChange->ProposeLocalEnergyDeposit(photonEnergy0 - photonEnergy1);
    return;
  }
  G4double eTotalEnergy = eKineticEnergy + electron_mass_c2;

  G4double electronE = photonEnergy0 * (1. - epsilon) + electron_mass_c2; 
  G4double electronP2 = electronE*electronE - electron_mass_c2*electron_mass_c2;
  G4double sinThetaE = -1.;
  G4double cosThetaE = 0.;
  if (electronP2 > 0.)
    {
      cosThetaE = (eTotalEnergy + photonEnergy1 )* (1. - epsilon) / std::sqrt(electronP2);
      sinThetaE = -1. * sqrt(1. - cosThetaE * cosThetaE); 
    }
  
  G4double eDirX = sinThetaE * std::cos(phi);
  G4double eDirY = sinThetaE * std::sin(phi);
  G4double eDirZ = cosThetaE;

  G4ThreeVector eDirection(eDirX,eDirY,eDirZ);
  eDirection.rotateUz(photonDirection0);

  // SI - The range test has been removed wrt original G4LowEnergyCompton class

  fParticleChange->ProposeLocalEnergyDeposit(bindingE);
  
  G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(),
						 eDirection,eKineticEnergy) ;
  fvect->push_back(dp);
}

