//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4LowEnergyCompton.cc,v 1.33 2001-11-07 20:47:29 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
//
// -------------------------------------------------------------------

#include "G4LowEnergyCompton.hh"
#include "Randomize.hh"
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
#include "G4VCrossSectionHandler.hh"
#include "G4CrossSectionHandler.hh"
#include "G4VEMDataSet.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"
#include "G4VRangeTest.hh"
#include "G4RangeTest.hh"

#include "G4CutsPerMaterialWarning.hh"

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
      G4Exception("G4LowEnergyCompton::G4LowEnergyCompton - energy outside intrinsic process validity range");
    }

  crossSectionHandler = new G4CrossSectionHandler;

  G4VDataSetAlgorithm* scatterInterpolation = new G4LogLogInterpolation;
  G4String scatterFile = "comp/ce-sf-";
  scatterFunctionData = new G4CompositeEMDataSet(scatterFile,scatterInterpolation,1.,1.);

  meanFreePathTable = 0;

  rangeTest = new G4RangeTest;

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

void G4LowEnergyCompton::BuildPhysicsTable(const G4ParticleDefinition& photon)
{
  
  G4CutsPerMaterialWarning warning;
  warning.PrintWarning(&photon);

  crossSectionHandler->Clear();
  G4String crossSectionFile = "comp/ce-cs-";
  crossSectionHandler->LoadData(crossSectionFile);

  delete meanFreePathTable;
  meanFreePathTable = crossSectionHandler->BuildMeanFreePathForMaterials();
}

G4VParticleChange* G4LowEnergyCompton::PostStepDoIt(const G4Track& aTrack, 
						    const G4Step&  aStep)
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
      aParticleChange.SetStatusChange(fStopAndKill);
      aParticleChange.SetEnergyChange(0.);
      aParticleChange.SetLocalEnergyDeposit(photonEnergy0);
      return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
    }

  G4double e0m = photonEnergy0 / electron_mass_c2 ;
  G4ParticleMomentum photonDirection0 = incidentPhoton->GetMomentumDirection();

  // Select randomly one element in the current material
  G4Material* material = aTrack.GetMaterial();
  G4int Z = crossSectionHandler->SelectRandomAtom(material,photonEnergy0);

  G4double epsilon0 = 1. / (1. + 2. * e0m);
  G4double epsilon0Sq = epsilon0 * epsilon0;
  G4double alpha1 = -log(epsilon0);
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
	  epsilon = exp(-alpha1 * G4UniformRand());  // pow(epsilon0,G4UniformRand())
	  epsilonSq = epsilon * epsilon; 
	}
      else
	{
	  epsilonSq = epsilon0Sq + (1. - epsilon0Sq) * G4UniformRand();
	  epsilon = sqrt(epsilonSq);
	}
      
      oneCosT = (1. - epsilon) / ( epsilon * e0m);
      sinT2 = oneCosT * (2. - oneCosT);      
      G4double x = sqrt(oneCosT/2.) / (wlPhoton/cm);
      G4double scatteringFunction = scatterFunctionData->FindValue(x,Z-1);
      gReject = (1. - epsilon * sinT2 / (1. + epsilonSq)) * scatteringFunction;
    
    }  while(gReject < G4UniformRand()*Z);

  G4double cosTheta = 1. - oneCosT;
  G4double sinTheta = sqrt (sinT2);
  G4double phi = twopi * G4UniformRand() ;
  G4double dirx = sinTheta * cos(phi);
  G4double diry = sinTheta * sin(phi);
  G4double dirz = cosTheta ;

  // Update G4VParticleChange for the scattered photon 
  
  G4ThreeVector photonDirection1(dirx,diry,dirz);
  photonDirection1.rotateUz(photonDirection0);
  aParticleChange.SetMomentumChange(photonDirection1) ;
  G4double photonEnergy1 = epsilon * photonEnergy0;

  if (photonEnergy1 > 0.)
    {
      aParticleChange.SetEnergyChange(photonEnergy1) ;
    }
  else
    {    
      aParticleChange.SetEnergyChange(0.) ;
      aParticleChange.SetStatusChange(fStopAndKill);
    }
  
  // Kinematics of the scattered electron 
  G4double eKineticEnergy = photonEnergy0 - photonEnergy1;

  // Generate the electron only if with large enough range w.r.t. cuts and safety

  G4double safety = aStep.GetPostStepPoint()->GetSafety();

  if (rangeTest->Escape(G4Electron::Electron(),material,eKineticEnergy,safety))
    {
      G4double eMomentum = sqrt(eKineticEnergy*(eKineticEnergy+2.*electron_mass_c2));
      G4ThreeVector eDirection((photonEnergy0 * photonDirection0 - 
				photonEnergy1 * photonDirection1) * (1./eMomentum));  
      G4DynamicParticle* electron = new G4DynamicParticle (G4Electron::Electron(),
							   eDirection,eKineticEnergy) ;
      aParticleChange.SetNumberOfSecondaries(1);
      aParticleChange.AddSecondary(electron);
      aParticleChange.SetLocalEnergyDeposit(0.); 
    }
  else
    {
      aParticleChange.SetNumberOfSecondaries(0);
      aParticleChange.SetLocalEnergyDeposit(eKineticEnergy);
    }

  return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep);
}

G4bool G4LowEnergyCompton::IsApplicable(const G4ParticleDefinition& particle)
{
  return ( &particle == G4Gamma::Gamma() ); 
}

G4double G4LowEnergyCompton::GetMeanFreePath(const G4Track& track, 
					     G4double previousStepSize, 
					     G4ForceCondition*)
{
  const G4DynamicParticle* photon = track.GetDynamicParticle();
  G4double energy = photon->GetKineticEnergy();
  G4Material* material = track.GetMaterial();
  size_t materialIndex = material->GetIndex();

  G4double meanFreePath;
  if (energy > highEnergyLimit) meanFreePath = meanFreePathTable->FindValue(highEnergyLimit,materialIndex);
  else if (energy < lowEnergyLimit) meanFreePath = DBL_MAX;
  else meanFreePath = meanFreePathTable->FindValue(energy,materialIndex);
  return meanFreePath;
}
