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
// --------------------------------------------------------------------
//
// $Id: G4LowEnergyRayleigh.cc,v 1.24 2001-08-28 16:05:20 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: A. Forti
//         Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -------- 
// Added Livermore data table construction methods A. Forti
// Added BuildMeanFreePath A. Forti
// Added PostStepDoIt A. Forti
// Added SelectRandomAtom A. Forti
// Added map of the elements  A.Forti
// 24.04.01 V.Ivanchenko remove RogueWave 
// 11.08.2001 MGP - Major revision according to a design iteration
//
// --------------------------------------------------------------------

#include "G4LowEnergyRayleigh.hh"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ForceCondition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4DynamicParticle.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleMomentum.hh"
#include "G4ThreeVector.hh"
#include "G4EnergyLossTables.hh"
#include "G4CrossSectionHandler.hh"
#include "G4VEMDataSet.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"


G4LowEnergyRayleigh::G4LowEnergyRayleigh(const G4String& processName)
  : G4VDiscreteProcess(processName),
    lowEnergyLimit(250*eV),             
    highEnergyLimit(100*GeV),
    intrinsicLowEnergyLimit(10*eV),
    intrinsicHighEnergyLimit(100*GeV)
{
  if (lowEnergyLimit < intrinsicLowEnergyLimit || 
      highEnergyLimit > intrinsicHighEnergyLimit)
    {
      G4Exception("G4LowEnergyRayleigh::G4LowEnergyRayleigh - energy limit outside intrinsic process validity range");
    }
  
  // The following pointer is owned by G4DataHandler
  G4VDataSetAlgorithm* crossSectionInterpolation = new G4LogLogInterpolation;
  crossSectionHandler = new G4CrossSectionHandler(crossSectionInterpolation);

  // The following pointer is owned by the process
  ffInterpolation = new G4LogLogInterpolation;
  G4String formFactorFile = "rayl/re-ff-";
  formFactorData = new G4CompositeEMDataSet(formFactorFile,ffInterpolation,1.,1.);

  meanFreePathTable = 0;

   if (verboseLevel > 0) 
     {
       G4cout << GetProcessName() << " is created " << G4endl
	      << "Energy range: " 
	      << lowEnergyLimit / keV << " keV - "
	      << highEnergyLimit / GeV << " GeV" 
	      << G4endl;
     }
}
 
G4LowEnergyRayleigh::~G4LowEnergyRayleigh()
{
  delete meanFreePathTable;
  delete crossSectionHandler;
  delete formFactorData;
  delete ffInterpolation;
}

void G4LowEnergyRayleigh::BuildPhysicsTable(const G4ParticleDefinition& photon)
{
  crossSectionHandler->Clear();
  G4String crossSectionFile = "rayl/re-cs-";
  crossSectionHandler->LoadData(crossSectionFile);

  delete meanFreePathTable;
  meanFreePathTable = crossSectionHandler->BuildMeanFreePathForMaterials();
}

G4VParticleChange* G4LowEnergyRayleigh::PostStepDoIt(const G4Track& aTrack, 
						     const G4Step& aStep)
{
  // The scattered gamma energy is sampled according to Form Factors 
  // multiplied by the Rayleigh distribution with a pure rejection technique.  
  // EGS4 W.R. Nelson et al. The EGS4 Code System. SLAC-Report-265 , December 1985 
  // Expression of the angular distribution as Rayleigh distribution and 
  // Form factors is taken from D. E. Cullen "A simple model of photon transport" 
  // NIM B Phys. 101 (1995). Method of sampling with form factors is different.
  // Reference to the article is from J. Stepanek New Photon, Positron
  // and Electron Interaction Data for GEANT in Energy Range from 1 eV to 10 TeV 
  // (draft). 

  aParticleChange.Initialize(aTrack);

  const G4DynamicParticle* incidentPhoton = aTrack.GetDynamicParticle();
  G4double photonEnergy0 = incidentPhoton->GetKineticEnergy();
  
  if (photonEnergy0 <= lowEnergyLimit)
    {
      aParticleChange.SetStatusChange(fStopAndKill);
      aParticleChange.SetEnergyChange(0.);
      aParticleChange.SetLocalEnergyDeposit(photonEnergy0);
      return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
    }

  //  G4double e0m = photonEnergy0 / electron_mass_c2 ;
  G4ParticleMomentum photonDirection0 = incidentPhoton->GetMomentumDirection();

  // Select randomly one element in the current material
  G4Material* material = aTrack.GetMaterial();
  G4int Z = crossSectionHandler->SelectRandomAtom(material,photonEnergy0);
  
  // Sample the energy of the scattered photon 

  G4double wlPhoton = h_Planck*c_light/photonEnergy0;

  G4double gReject;
  G4double randomFormFactor;
  G4double cosTheta;
  G4double sinTheta;

  do
    {
      G4double thetaHalf = G4UniformRand() * pi / 2.;
      G4double sinThetaHalf = sin(thetaHalf);
      G4double x = sinThetaHalf / (wlPhoton/cm);
      G4double dataFormFactor = formFactorData->FindValue(x,Z-1);
      randomFormFactor = G4UniformRand() * Z * Z;
      G4double theta = thetaHalf*2;
      cosTheta = cos(theta);
      sinTheta = sin(theta);
      G4double sqrRayl = 1 + cosTheta * cosTheta;    
      gReject = sqrRayl * dataFormFactor * dataFormFactor;

    } while( gReject < randomFormFactor);
  
  // Scattered photon angles. ( Z - axis along the parent photon)
  G4double phi = twopi * G4UniformRand() ;
  G4double dirX = sinTheta*cos(phi);
  G4double dirY = sinTheta*sin(phi);
  G4double dirZ = cosTheta;
  
  // Update G4VParticleChange for the scattered photon 
  G4ThreeVector photonDirection1(dirX, dirY, dirZ);

  photonDirection1.rotateUz(photonDirection0);
  aParticleChange.SetEnergyChange(photonEnergy0);
  aParticleChange.SetMomentumChange(photonDirection1);
  
  aParticleChange.SetNumberOfSecondaries(0);

#ifdef G4VERBOSE
  if(verboseLevel > 15) G4cout << "LE Rayleigh PostStepDoIt" << G4endl;
#endif

  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

G4bool G4LowEnergyRayleigh::IsApplicable(const G4ParticleDefinition& particle)
{
  return ( &particle == G4Gamma::Gamma() ); 
}

G4double G4LowEnergyRayleigh::GetMeanFreePath(const G4Track& track, 
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
