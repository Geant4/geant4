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
// --------------------------------------------------------------------
///
//
// 
// --------------------------------------------------------------
//
// Author: A. Forti
//         Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -------- 
// 02/03/1999 A. Forti 1st implementation
// 14.03.2000 Veronique Lefebure;
// Change initialisation of lowestEnergyLimit from 1.22 to 1.022.
// Note that the hard coded value 1.022 should be used instead of
// 2*electron_mass_c2 in order to agree with the value of the data bank EPDL97
// 24.04.01 V.Ivanchenko remove RogueWave
// 27.07.01 F.Longo correct bug in energy distribution
// 21.01.03 V.Ivanchenko Cut per region
// 25.03.03 F.Longo fix in angular distribution of e+/e-
// 24.04.03 V.Ivanchenko - Cut per region mfpt
//
// --------------------------------------------------------------

#include "G4LowEnergyGammaConversion.hh"

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
#include "G4Positron.hh"
#include "G4IonisParamElm.hh"
#include "G4Material.hh"
#include "G4RDVCrossSectionHandler.hh"
#include "G4RDCrossSectionHandler.hh"
#include "G4RDVEMDataSet.hh"
#include "G4RDVDataSetAlgorithm.hh"
#include "G4RDLogLogInterpolation.hh"
#include "G4RDVRangeTest.hh"
#include "G4RDRangeTest.hh"
#include "G4MaterialCutsCouple.hh"

G4LowEnergyGammaConversion::G4LowEnergyGammaConversion(const G4String& processName)
  : G4VDiscreteProcess(processName),
    lowEnergyLimit(1.022000*MeV),
    highEnergyLimit(100*GeV),
    intrinsicLowEnergyLimit(1.022000*MeV),
    intrinsicHighEnergyLimit(100*GeV),
    smallEnergy(2.*MeV)

{
  if (lowEnergyLimit < intrinsicLowEnergyLimit || 
      highEnergyLimit > intrinsicHighEnergyLimit)
    {
      G4Exception("G4LowEnergyGammaConversion::G4LowEnergyGammaConversion()",
                  "OutOfRange", FatalException,
                  "Energy limit outside intrinsic process validity range!");
    }

  // The following pointer is owned by G4DataHandler
  
  crossSectionHandler = new G4RDCrossSectionHandler();
  crossSectionHandler->Initialise(0,1.0220*MeV,100.*GeV,400);
  meanFreePathTable = 0;
  rangeTest = new G4RDRangeTest;

   if (verboseLevel > 0) 
     {
       G4cout << GetProcessName() << " is created " << G4endl
	      << "Energy range: " 
	      << lowEnergyLimit / MeV << " MeV - "
	      << highEnergyLimit / GeV << " GeV" 
	      << G4endl;
     }
}
 
G4LowEnergyGammaConversion::~G4LowEnergyGammaConversion()
{
  delete meanFreePathTable;
  delete crossSectionHandler;
  delete rangeTest;
}

void G4LowEnergyGammaConversion::BuildPhysicsTable(const G4ParticleDefinition& )
{

  crossSectionHandler->Clear();
  G4String crossSectionFile = "pair/pp-cs-";
  crossSectionHandler->LoadData(crossSectionFile);

  delete meanFreePathTable;
  meanFreePathTable = crossSectionHandler->BuildMeanFreePathForMaterials();
}

G4VParticleChange* G4LowEnergyGammaConversion::PostStepDoIt(const G4Track& aTrack,
							    const G4Step& aStep)
{
// The energies of the e+ e- secondaries are sampled using the Bethe - Heitler
// cross sections with Coulomb correction. A modified version of the random
// number techniques of Butcher & Messel is used (Nuc Phys 20(1960),15).

// Note 1 : Effects due to the breakdown of the Born approximation at low
// energy are ignored.
// Note 2 : The differential cross section implicitly takes account of
// pair creation in both nuclear and atomic electron fields. However triplet
// prodution is not generated.

  aParticleChange.Initialize(aTrack);

  const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();

  const G4DynamicParticle* incidentPhoton = aTrack.GetDynamicParticle();
  G4double photonEnergy = incidentPhoton->GetKineticEnergy();
  G4ParticleMomentum photonDirection = incidentPhoton->GetMomentumDirection();

  G4double epsilon ;
  G4double epsilon0 = electron_mass_c2 / photonEnergy ;

  // Do it fast if photon energy < 2. MeV
  if (photonEnergy < smallEnergy )
    {
      epsilon = epsilon0 + (0.5 - epsilon0) * G4UniformRand();
    }
  else
    {
      // Select randomly one element in the current material
      const G4Element* element = crossSectionHandler->SelectRandomElement(couple,photonEnergy);

      if (element == 0)
	{
	  G4cout << "G4LowEnergyGammaConversion::PostStepDoIt - element = 0" << G4endl;
	}
      G4IonisParamElm* ionisation = element->GetIonisation();
       if (ionisation == 0)
	{
	  G4cout << "G4LowEnergyGammaConversion::PostStepDoIt - ionisation = 0" << G4endl;
	}

      // Extract Coulomb factor for this Element
      G4double fZ = 8. * (ionisation->GetlogZ3());
      if (photonEnergy > 50. * MeV) fZ += 8. * (element->GetfCoulomb());

      // Limits of the screening variable
      G4double screenFactor = 136. * epsilon0 / (element->GetIonisation()->GetZ3()) ;
      G4double screenMax = std::exp ((42.24 - fZ)/8.368) - 0.952 ;
      G4double screenMin = std::min(4.*screenFactor,screenMax) ;

      // Limits of the energy sampling
      G4double epsilon1 = 0.5 - 0.5 * std::sqrt(1. - screenMin / screenMax) ;
      G4double epsilonMin = std::max(epsilon0,epsilon1);
      G4double epsilonRange = 0.5 - epsilonMin ;

      // Sample the energy rate of the created electron (or positron)
      G4double screen;
      G4double gReject ;

      G4double f10 = ScreenFunction1(screenMin) - fZ;
      G4double f20 = ScreenFunction2(screenMin) - fZ;
      G4double normF1 = std::max(f10 * epsilonRange * epsilonRange,0.);
      G4double normF2 = std::max(1.5 * f20,0.);

      do {
	if (normF1 / (normF1 + normF2) > G4UniformRand() )
	  {
	    epsilon = 0.5 - epsilonRange * std::pow(G4UniformRand(), 0.3333) ;
	    screen = screenFactor / (epsilon * (1. - epsilon));
	    gReject = (ScreenFunction1(screen) - fZ) / f10 ;
	  }
	else
	  {
	    epsilon = epsilonMin + epsilonRange * G4UniformRand();
	    screen = screenFactor / (epsilon * (1 - epsilon));
	    gReject = (ScreenFunction2(screen) - fZ) / f20 ;
	  }
      } while ( gReject < G4UniformRand() );

    }   //  End of epsilon sampling

  // Fix charges randomly

  G4double electronTotEnergy;
  G4double positronTotEnergy;

  if (CLHEP::RandBit::shootBit())
    {
      electronTotEnergy = (1. - epsilon) * photonEnergy;
      positronTotEnergy = epsilon * photonEnergy;
    }
  else
    {
      positronTotEnergy = (1. - epsilon) * photonEnergy;
      electronTotEnergy = epsilon * photonEnergy;
    }

  // Scattered electron (positron) angles. ( Z - axis along the parent photon)
  // Universal distribution suggested by L. Urban (Geant3 manual (1993) Phys211),
  // derived from Tsai distribution (Rev. Mod. Phys. 49, 421 (1977)

  G4double u;
  const G4double a1 = 0.625;
  G4double a2 = 3. * a1;
  //  G4double d = 27. ;

  //  if (9. / (9. + d) > G4UniformRand())
  if (0.25 > G4UniformRand())
    {
      u = - std::log(G4UniformRand() * G4UniformRand()) / a1 ;
    }
  else
    {
      u = - std::log(G4UniformRand() * G4UniformRand()) / a2 ;
    }

  G4double thetaEle = u*electron_mass_c2/electronTotEnergy;
  G4double thetaPos = u*electron_mass_c2/positronTotEnergy;
  G4double phi  = twopi * G4UniformRand();

  G4double dxEle= std::sin(thetaEle)*std::cos(phi),dyEle= std::sin(thetaEle)*std::sin(phi),dzEle=std::cos(thetaEle);
  G4double dxPos=-std::sin(thetaPos)*std::cos(phi),dyPos=-std::sin(thetaPos)*std::sin(phi),dzPos=std::cos(thetaPos);
  
  
  // Kinematics of the created pair:
  // the electron and positron are assumed to have a symetric angular 
  // distribution with respect to the Z axis along the parent photon
  
  G4double localEnergyDeposit = 0. ;
  
  aParticleChange.SetNumberOfSecondaries(2) ; 
  G4double electronKineEnergy = std::max(0.,electronTotEnergy - electron_mass_c2) ;
  
  // Generate the electron only if with large enough range w.r.t. cuts and safety
  
  G4double safety = aStep.GetPostStepPoint()->GetSafety();
  
  if (rangeTest->Escape(G4Electron::Electron(),couple,electronKineEnergy,safety))
    {
      G4ThreeVector electronDirection (dxEle, dyEle, dzEle);
      electronDirection.rotateUz(photonDirection);
      
      G4DynamicParticle* particle1 = new G4DynamicParticle (G4Electron::Electron(),
							    electronDirection,
							    electronKineEnergy);
      aParticleChange.AddSecondary(particle1) ;
    }
  else
    {
      localEnergyDeposit += electronKineEnergy ;
    }

  // The e+ is always created (even with kinetic energy = 0) for further annihilation
  G4double positronKineEnergy = std::max(0.,positronTotEnergy - electron_mass_c2) ;

  // Is the local energy deposit correct, if the positron is always created?
  if (! (rangeTest->Escape(G4Positron::Positron(),couple,positronKineEnergy,safety)))
    {
      localEnergyDeposit += positronKineEnergy ;
      positronKineEnergy = 0. ;
    }

  G4ThreeVector positronDirection (dxPos, dyPos, dzPos);
  positronDirection.rotateUz(photonDirection);   
  
  // Create G4DynamicParticle object for the particle2 
  G4DynamicParticle* particle2 = new G4DynamicParticle(G4Positron::Positron(),
						       positronDirection, positronKineEnergy);
  aParticleChange.AddSecondary(particle2) ; 

  aParticleChange.ProposeLocalEnergyDeposit(localEnergyDeposit) ;
  
  // Kill the incident photon 
  aParticleChange.ProposeMomentumDirection(0.,0.,0.) ;
  aParticleChange.ProposeEnergy(0.) ; 
  aParticleChange.ProposeTrackStatus(fStopAndKill) ;

  //  Reset NbOfInteractionLengthLeft and return aParticleChange
  return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
}

G4bool G4LowEnergyGammaConversion::IsApplicable(const G4ParticleDefinition& particle)
{
  return ( &particle == G4Gamma::Gamma() ); 
}

G4double G4LowEnergyGammaConversion::GetMeanFreePath(const G4Track& track, 
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

G4double G4LowEnergyGammaConversion::ScreenFunction1(G4double screenVariable)
{
  // Compute the value of the screening function 3*phi1 - phi2

  G4double value;
  
  if (screenVariable > 1.)
    value = 42.24 - 8.368 * std::log(screenVariable + 0.952);
  else
    value = 42.392 - screenVariable * (7.796 - 1.961 * screenVariable);
  
  return value;
} 

G4double G4LowEnergyGammaConversion::ScreenFunction2(G4double screenVariable)
{
  // Compute the value of the screening function 1.5*phi1 - 0.5*phi2
  
  G4double value;
  
  if (screenVariable > 1.)
    value = 42.24 - 8.368 * std::log(screenVariable + 0.952);
  else
    value = 41.405 - screenVariable * (5.828 - 0.8945 * screenVariable);
  
  return value;
} 
