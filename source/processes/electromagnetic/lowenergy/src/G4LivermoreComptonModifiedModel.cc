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
// $Id: G4LivermoreComptonModifiedModel.cc 95950 2016-03-03 10:42:48Z gcosmo $
//
//
// Author: Sebastien Incerti
//         30 October 2008
//         on base of G4LowEnergyCompton developed by A.Forti and M.G.Pia
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
// 30 May 2011   V Ivanchenko Migration to model design for deexcitation

#include "G4LivermoreComptonModifiedModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4LossTableManager.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4AtomicShell.hh"
#include "G4CrossSectionHandler.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4LogLogInterpolation.hh"
#include "G4Gamma.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermoreComptonModifiedModel::G4LivermoreComptonModifiedModel(const G4ParticleDefinition*,
						 const G4String& nam)
  :G4VEmModel(nam),fParticleChange(0),isInitialised(false),
   scatterFunctionData(0),
   crossSectionHandler(0),fAtomDeexcitation(0)
{
  verboseLevel=0 ;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if(  verboseLevel>0 )
    G4cout << "Livermore Modified Compton model is constructed " << G4endl;

  //Mark this model as "applicable" for atomic deexcitation
  SetDeexcitationFlag(true);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermoreComptonModifiedModel::~G4LivermoreComptonModifiedModel()
{
  delete crossSectionHandler;
  delete scatterFunctionData;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreComptonModifiedModel::Initialise(const G4ParticleDefinition* particle,
					 const G4DataVector& cuts)
{
  if (verboseLevel > 2) {
    G4cout << "Calling G4LivermoreComptonModifiedModel::Initialise()" << G4endl;
  }

  if (crossSectionHandler)
  {
    crossSectionHandler->Clear();
    delete crossSectionHandler;
  }
  delete scatterFunctionData;

  // Reading of data files - all materials are read
  crossSectionHandler = new G4CrossSectionHandler;
  G4String crossSectionFile = "comp/ce-cs-";
  crossSectionHandler->LoadData(crossSectionFile);

  G4VDataSetAlgorithm* scatterInterpolation = new G4LogLogInterpolation;
  G4String scatterFile = "comp/ce-sf-";
  scatterFunctionData = new G4CompositeEMDataSet(scatterInterpolation, 1., 1.);
  scatterFunctionData->LoadData(scatterFile);

  // For Doppler broadening
  shellData.SetOccupancyData();
  G4String file = "/doppler/shell-doppler";
  shellData.LoadData(file);

  InitialiseElementSelectors(particle,cuts);

  if (verboseLevel > 2) {
    G4cout << "Loaded cross section files for Livermore Modified Compton model" << G4endl;
  }

  if(isInitialised) { return; }
  isInitialised = true;

  fParticleChange = GetParticleChangeForGamma();

  fAtomDeexcitation  = G4LossTableManager::Instance()->AtomDeexcitation();

  if(  verboseLevel>0 ) {
    G4cout << "Livermore modified Compton model is initialized " << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / eV << " eV - "
	   << HighEnergyLimit() / GeV << " GeV"
	   << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermoreComptonModifiedModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double GammaEnergy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
  if (verboseLevel > 3) {
    G4cout << "Calling ComputeCrossSectionPerAtom() of G4LivermoreComptonModifiedModel" << G4endl;
  }
  if (GammaEnergy < LowEnergyLimit())
    { return 0.0; }

  G4double cs = crossSectionHandler->FindValue(G4int(Z), GammaEnergy);
  return cs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreComptonModifiedModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
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
    G4cout << "G4LivermoreComptonModifiedModel::SampleSecondaries() E(MeV)= "
	   << photonEnergy0/MeV << " in " << couple->GetMaterial()->GetName()
	   << G4endl;
  }

  // do nothing below the threshold
  // should never get here because the XS is zero below the limit
  if (photonEnergy0 < LowEnergyLimit())
    return ;

  G4double e0m = photonEnergy0 / electron_mass_c2 ;
  G4ParticleMomentum photonDirection0 = aDynamicGamma->GetMomentumDirection();

  // Select randomly one element in the current material
  const G4ParticleDefinition* particle =  aDynamicGamma->GetDefinition();
  const G4Element* elm = SelectRandomAtom(couple,particle,photonEnergy0);
  G4int Z = (G4int)elm->GetZ();

  G4double epsilon0Local = 1. / (1. + 2. * e0m);
  G4double epsilon0Sq = epsilon0Local * epsilon0Local;
  G4double alpha1 = -std::log(epsilon0Local);
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
	  // std::pow(epsilon0Local,G4UniformRand())
	  epsilon = G4Exp(-alpha1 * G4UniformRand());
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
  G4double systemE = 0;
  G4double ePAU = -1;
  G4int shellIdx = 0;
  G4double vel_c = 299792458;
  G4double momentum_au_to_nat = 1.992851740*std::pow(10.,-24.);
  G4double e_mass_kg = 9.10938188 * std::pow(10.,-31.);
  G4double eMax = -1;
  G4double Alpha=0;
  do
    {
      ++iteration;
      // Select shell based on shell occupancy
      shellIdx = shellData.SelectRandomShell(Z);
      bindingE = shellData.BindingEnergy(Z,shellIdx);



      // Randomly sample bound electron momentum
      // (memento: the data set is in Atomic Units)
      G4double pSample = profileData.RandomSelectMomentum(Z,shellIdx);
      // Rescale from atomic units


      //Kinetic energy of target electron


      // Reverse vector projection onto scattering vector

      do {
         Alpha = G4UniformRand()*pi/2.0;
         } while(Alpha >= (pi/2.0));

      ePAU = pSample / std::cos(Alpha);

      // Convert to SI and the calculate electron energy in natural units

      G4double ePSI = ePAU * momentum_au_to_nat;
      G4double u_temp = sqrt( ((ePSI*ePSI)*(vel_c*vel_c)) / ((e_mass_kg*e_mass_kg)*(vel_c*vel_c)+(ePSI*ePSI)))/vel_c;
      G4double eEIncident = electron_mass_c2 / sqrt( 1 - (u_temp*u_temp));

      //Total energy of the system
      systemE = eEIncident+photonEnergy0;

      eMax = systemE - bindingE - electron_mass_c2;
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
 	  if (G4UniformRand() < 0.5) { photonE = (var4 - varSqrt) * scale; }
	  else                       { photonE = (var4 + varSqrt) * scale; }
	}
      else
	{
	  photonE = -1.;
	}
    } while ( iteration <= maxDopplerIterations &&
	     (photonE < 0. || photonE > eMax ) );

  // End of recalculation of photon energy with Doppler broadening
  // Kinematics of the scattered electron
  G4double eKineticEnergy = systemE - photonE - bindingE - electron_mass_c2;

  // protection against negative final energy: no e- is created
   G4double eDirX = 0.0;
   G4double eDirY = 0.0;
   G4double eDirZ = 1.0;

  if(eKineticEnergy < 0.0) {
    G4cout << "Error, kinetic energy of electron less than zero" << G4endl;
    }

 else{
    // Estimation of Compton electron polar angle taken from:
    // The EGSnrc Code System: Monte Carlo Simulation of Electron and Photon Transport
    // Eqn 2.2.25 Pg 42, NRCC Report PIRS-701
    G4double E_num = photonEnergy0 - photonE*cosTheta;
    G4double E_dom = sqrt(photonEnergy0*photonEnergy0 + photonE*photonE -2*photonEnergy0*photonE*cosTheta);
    G4double cosThetaE = E_num / E_dom;
    G4double sinThetaE = -sqrt((1. - cosThetaE) * (1. + cosThetaE));

    eDirX = sinThetaE * std::cos(phi);
    eDirY = sinThetaE * std::sin(phi);
    eDirZ = cosThetaE;

    G4ThreeVector eDirection(eDirX,eDirY,eDirZ);
    eDirection.rotateUz(photonDirection0);
    G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(),
                                                   eDirection,eKineticEnergy) ;
    fvect->push_back(dp);
   }


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

        if (iteration < maxDopplerIterations)
        {
         G4ThreeVector eDirection(eDirX,eDirY,eDirZ);
         eDirection.rotateUz(photonDirection0);
         G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(),
                                                       eDirection,eKineticEnergy) ;
         fvect->push_back(dp);
        }
    }
  else
    {
      photonEnergy1 = 0.;
      fParticleChange->SetProposedKineticEnergy(0.) ;
      fParticleChange->ProposeTrackStatus(fStopAndKill);
     }

  // sample deexcitation
  //
  if(fAtomDeexcitation && iteration < maxDopplerIterations) {
    G4int index = couple->GetIndex();
    if(fAtomDeexcitation->CheckDeexcitationActiveRegion(index)) {
      size_t nbefore = fvect->size();
      G4AtomicShellEnumerator as = G4AtomicShellEnumerator(shellIdx);
      const G4AtomicShell* shell = fAtomDeexcitation->GetAtomicShell(Z, as);
      fAtomDeexcitation->GenerateParticles(fvect, shell, Z, index);
      size_t nafter = fvect->size();
      if(nafter > nbefore) {
	for (size_t i=nbefore; i<nafter; ++i) {
	  bindingE -= ((*fvect)[i])->GetKineticEnergy();
	}
      }
    }
  }
  if(bindingE < 0.0) { bindingE = 0.0; }
  fParticleChange->ProposeLocalEnergyDeposit(bindingE);
}
