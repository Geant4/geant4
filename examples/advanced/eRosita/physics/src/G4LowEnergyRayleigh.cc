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
//
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
// 06.10.2001 MGP - Added strategy to test range for secondary generation
// 05.06.2002 F.Longo and G.Depaola  - bug fixed in angular distribution
// 20.10.2002 G. Depaola   - Change sampling method of theta
// 22.01.2003 V.Ivanchenko - Cut per region
// 24.04.2003 V.Ivanchenko - Cut per region mfpt
//
// --------------------------------------------------------------------

#include "G4LowEnergyRayleigh.hh"
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
#include "G4RDVCrossSectionHandler.hh"
#include "G4RDCrossSectionHandler.hh"
#include "G4RDVEMDataSet.hh"
#include "G4RDCompositeEMDataSet.hh"
#include "G4RDVDataSetAlgorithm.hh"
#include "G4RDLogLogInterpolation.hh"

#include "G4MaterialCutsCouple.hh"

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
      G4Exception("G4LowEnergyRayleigh::G4LowEnergyRayleigh()",
                  "OutOfRange", FatalException,
                  "Energy limit outside intrinsic process validity range!");
    }

  crossSectionHandler = new G4RDCrossSectionHandler();

  G4RDVDataSetAlgorithm* ffInterpolation = new G4RDLogLogInterpolation;
  G4String formFactorFile = "rayl/re-ff-";
  formFactorData = new G4RDCompositeEMDataSet(ffInterpolation,1.,1.);
  formFactorData->LoadData(formFactorFile);

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
}

void G4LowEnergyRayleigh::BuildPhysicsTable(const G4ParticleDefinition& )
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

  aParticleChange.Initialize(aTrack);

  const G4DynamicParticle* incidentPhoton = aTrack.GetDynamicParticle();
  G4double photonEnergy0 = incidentPhoton->GetKineticEnergy();

  if (photonEnergy0 <= lowEnergyLimit)
    {
      aParticleChange.ProposeTrackStatus(fStopAndKill);
      aParticleChange.ProposeEnergy(0.);
      aParticleChange.ProposeLocalEnergyDeposit(photonEnergy0);
      return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
    }

  //  G4double e0m = photonEnergy0 / electron_mass_c2 ;
  G4ParticleMomentum photonDirection0 = incidentPhoton->GetMomentumDirection();

  // Select randomly one element in the current material
  const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
  G4int Z = crossSectionHandler->SelectRandomAtom(couple,photonEnergy0);

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

      G4double sinThetaHalf = std::sqrt((1. - cosTheta) / 2.);
      x = sinThetaHalf / (wlPhoton/cm);
      if (x > 1.e+005)
         dataFormFactor = formFactorData->FindValue(x,Z-1);
      else
         dataFormFactor = formFactorData->FindValue(0.,Z-1);
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
  aParticleChange.ProposeEnergy(photonEnergy0);
  aParticleChange.ProposeMomentumDirection(photonDirection1);

  aParticleChange.SetNumberOfSecondaries(0);

  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

G4bool G4LowEnergyRayleigh::IsApplicable(const G4ParticleDefinition& particle)
{
  return ( &particle == G4Gamma::Gamma() ); 
}

G4double G4LowEnergyRayleigh::GetMeanFreePath(const G4Track& track, 
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
