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
// Author: A. Forti
//         Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// --------
// Added Livermore data table construction methods A. Forti
// Modified BuildMeanFreePath to read new data tables A. Forti
// Modified PostStepDoIt to insert sampling with EPDL97 data A. Forti
// Added SelectRandomAtom A. Forti
// Added map of the elements A. Forti
// 24.04.2001 V.Ivanchenko - Remove RogueWave
// 06.08.2001 MGP          - Revised according to a design iteration
// 22.01.2003 V.Ivanchenko - Cut per region
// 10.03.2003 V.Ivanchenko - Remove CutPerMaterial warning
// 24.04.2003 V.Ivanchenko - Cut per region mfpt
//
// -------------------------------------------------------------------

#include "G4LowEnergyCompton.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ForceCondition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4DynamicParticle.hh"
#include "G4VParticleChange.hh"
#include "G4ThreeVector.hh"
#include "G4EnergyLossTables.hh"
#include "G4RDVCrossSectionHandler.hh"
#include "G4RDCrossSectionHandler.hh"
#include "G4RDVEMDataSet.hh"
#include "G4RDCompositeEMDataSet.hh"
#include "G4RDVDataSetAlgorithm.hh"
#include "G4RDLogLogInterpolation.hh"
#include "G4RDVRangeTest.hh"
#include "G4RDRangeTest.hh"
#include "G4RDRangeNoTest.hh"
#include "G4MaterialCutsCouple.hh"


G4LowEnergyCompton::G4LowEnergyCompton(const G4String& processName)
  : G4VDiscreteProcess(processName),
    lowEnergyLimit(250*eV),
    highEnergyLimit(100*GeV),
    intrinsicLowEnergyLimit(10*eV),
    intrinsicHighEnergyLimit(100*GeV)
{
  if (lowEnergyLimit < intrinsicLowEnergyLimit ||
      highEnergyLimit > intrinsicHighEnergyLimit)
    {
      G4Exception("G4LowEnergyCompton::G4LowEnergyCompton()",
                  "OutOfRange", FatalException,
                  "Energy outside intrinsic process validity range!");
    }

  crossSectionHandler = new G4RDCrossSectionHandler;

  G4RDVDataSetAlgorithm* scatterInterpolation = new G4RDLogLogInterpolation;
  G4String scatterFile = "comp/ce-sf-";
  scatterFunctionData = new G4RDCompositeEMDataSet(scatterInterpolation, 1., 1.);
  scatterFunctionData->LoadData(scatterFile);

  meanFreePathTable = 0;

  rangeTest = new G4RDRangeNoTest;

  // For Doppler broadening
  shellData.SetOccupancyData();

   if (verboseLevel > 0)
     {
       G4cout << GetProcessName() << " is created " << G4endl
	      << "Energy range: "
	      << lowEnergyLimit / keV << " keV - "
	      << highEnergyLimit / GeV << " GeV"
	      << G4endl;
     }
}

G4LowEnergyCompton::~G4LowEnergyCompton()
{
  delete meanFreePathTable;
  delete crossSectionHandler;
  delete scatterFunctionData;
  delete rangeTest;
}

void G4LowEnergyCompton::BuildPhysicsTable(const G4ParticleDefinition& )
{

  crossSectionHandler->Clear();
  G4String crossSectionFile = "comp/ce-cs-";
  crossSectionHandler->LoadData(crossSectionFile);

  delete meanFreePathTable;
  meanFreePathTable = crossSectionHandler->BuildMeanFreePathForMaterials();

  // For Doppler broadening
  G4String file = "/doppler/shell-doppler";
  shellData.LoadData(file);
}

G4VParticleChange* G4LowEnergyCompton::PostStepDoIt(const G4Track& aTrack,
						    const G4Step& aStep)
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

  aParticleChange.Initialize(aTrack);

  // Dynamic particle quantities
  const G4DynamicParticle* incidentPhoton = aTrack.GetDynamicParticle();
  G4double photonEnergy0 = incidentPhoton->GetKineticEnergy();

  if (photonEnergy0 <= lowEnergyLimit)
    {
      aParticleChange.ProposeTrackStatus(fStopAndKill);
      aParticleChange.ProposeEnergy(0.);
      aParticleChange.ProposeLocalEnergyDeposit(photonEnergy0);
      return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
    }

  G4double e0m = photonEnergy0 / electron_mass_c2 ;
  G4ParticleMomentum photonDirection0 = incidentPhoton->GetMomentumDirection();
  G4double epsilon0 = 1. / (1. + 2. * e0m);
  G4double epsilon0Sq = epsilon0 * epsilon0;
  G4double alpha1 = -std::log(epsilon0);
  G4double alpha2 = 0.5 * (1. - epsilon0Sq);
  G4double wlPhoton = h_Planck*c_light/photonEnergy0;

  // Select randomly one element in the current material
  const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
  G4int Z = crossSectionHandler->SelectRandomAtom(couple,photonEnergy0);

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
	  epsilon = std::exp(-alpha1 * G4UniformRand());  // std::pow(epsilon0,G4UniformRand())
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

    }  while(gReject < G4UniformRand()*Z);

  G4double cosTheta = 1. - oneCosT;
  G4double sinTheta = std::sqrt (sinT2);
  G4double phi = twopi * G4UniformRand() ;
  G4double dirX = sinTheta * std::cos(phi);
  G4double dirY = sinTheta * std::sin(phi);
  G4double dirZ = cosTheta ;

  // Doppler broadening -  Method based on:
  // Y. Namito, S. Ban and H. Hirayama, 
  // "Implementation of the Doppler Broadening of a Compton-Scattered Photon Into the EGS4 Code" 
  // NIM A 349, pp. 489-494, 1994
  
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
     
      // Randomly sample bound electron momentum (memento: the data set is in Atomic Units)
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

  G4ThreeVector photonDirection1(dirX,dirY,dirZ);
  photonDirection1.rotateUz(photonDirection0);
  aParticleChange.ProposeMomentumDirection(photonDirection1);
 
  G4double photonEnergy1 = photonE;
  //G4cout << "--> PHOTONENERGY1 = " << photonE/keV << G4endl;

  if (photonEnergy1 > 0.)
    {
      aParticleChange.ProposeEnergy(photonEnergy1) ;
    }
  else
    {
      aParticleChange.ProposeEnergy(0.) ;
      aParticleChange.ProposeTrackStatus(fStopAndKill);
    }

  // Kinematics of the scattered electron
  G4double eKineticEnergy = photonEnergy0 - photonEnergy1 - bindingE;
  G4double eTotalEnergy = eKineticEnergy + electron_mass_c2;

  G4double electronE = photonEnergy0 * (1. - epsilon) + electron_mass_c2; 
  G4double electronP2 = electronE*electronE - electron_mass_c2*electron_mass_c2;
  G4double sinThetaE = -1.;
  G4double cosThetaE = 0.;
  if (electronP2 > 0.)
    {
      cosThetaE = (eTotalEnergy + photonEnergy1 )* (1. - epsilon) / std::sqrt(electronP2);
      sinThetaE = -1. * std::sqrt(1. - cosThetaE * cosThetaE); 
    }
  
  G4double eDirX = sinThetaE * std::cos(phi);
  G4double eDirY = sinThetaE * std::sin(phi);
  G4double eDirZ = cosThetaE;

  // Generate the electron only if with large enough range w.r.t. cuts and safety

  G4double safety = aStep.GetPostStepPoint()->GetSafety();

  if (rangeTest->Escape(G4Electron::Electron(),couple,eKineticEnergy,safety))
    {
      G4ThreeVector eDirection(eDirX,eDirY,eDirZ);
      eDirection.rotateUz(photonDirection0);

      G4DynamicParticle* electron = new G4DynamicParticle (G4Electron::Electron(),eDirection,eKineticEnergy) ;
      aParticleChange.SetNumberOfSecondaries(1);
      aParticleChange.AddSecondary(electron);
      // Binding energy deposited locally
      aParticleChange.ProposeLocalEnergyDeposit(bindingE);
    }
  else
    {
      aParticleChange.SetNumberOfSecondaries(0);
      aParticleChange.ProposeLocalEnergyDeposit(eKineticEnergy + bindingE);
    }

  return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep);
}

G4bool G4LowEnergyCompton::IsApplicable(const G4ParticleDefinition& particle)
{
  return ( &particle == G4Gamma::Gamma() );
}

G4double G4LowEnergyCompton::GetMeanFreePath(const G4Track& track,
					     G4double, // previousStepSize
					     G4ForceCondition*)
{
  const G4DynamicParticle* photon = track.GetDynamicParticle();
  G4double energy = photon->GetKineticEnergy();
  const G4MaterialCutsCouple* couple = track.GetMaterialCutsCouple();
  size_t materialIndex = couple->GetIndex();

  G4double meanFreePath;
  if (energy > highEnergyLimit) meanFreePath = meanFreePathTable->FindValue(highEnergyLimit,materialIndex);
  else if (energy < lowEnergyLimit) meanFreePath = DBL_MAX;
  else meanFreePath = meanFreePathTable->FindValue(energy,materialIndex);
  return meanFreePath;
}
