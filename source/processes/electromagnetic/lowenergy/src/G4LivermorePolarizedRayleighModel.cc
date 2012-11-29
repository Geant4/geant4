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
// Author: Sebastien Incerti
//         30 October 2008
//         on base of G4LowEnergyPolarizedRayleigh developed by R. Capra
//
// History:
// --------
// 02 May 2009   S Incerti as V. Ivanchenko proposed in G4LivermoreRayleighModel.cc
//
// Cleanup initialisation and generation of secondaries:
//                  - apply internal high-energy limit only in constructor 
//                  - do not apply low-energy limit (default is 0)
//                  - remove GetMeanFreePath method and table
//                  - remove initialisation of element selector 
//                  - use G4ElementSelector

#include "G4LivermorePolarizedRayleighModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePolarizedRayleighModel::G4LivermorePolarizedRayleighModel(const G4ParticleDefinition*,
                                             const G4String& nam)
  :G4VEmModel(nam),fParticleChange(0),isInitialised(false),
   crossSectionHandler(0),formFactorData(0)
{
  lowEnergyLimit = 250 * eV; 
  highEnergyLimit = 100 * GeV;
  
  //SetLowEnergyLimit(lowEnergyLimit);
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
    G4cout << "Livermore Polarized Rayleigh is constructed " << G4endl
         << "Energy range: "
         << lowEnergyLimit / eV << " eV - "
         << highEnergyLimit / GeV << " GeV"
         << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePolarizedRayleighModel::~G4LivermorePolarizedRayleighModel()
{  
  if (crossSectionHandler) delete crossSectionHandler;
  if (formFactorData) delete formFactorData;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedRayleighModel::Initialise(const G4ParticleDefinition* particle,
                                       const G4DataVector& cuts)
{
// Rayleigh process:                      The Quantum Theory of Radiation
//                                        W. Heitler,       Oxford at the Clarendon Press, Oxford (1954)                                                 
// Scattering function:                   A simple model of photon transport
//                                        D.E. Cullen,      Nucl. Instr. Meth. in Phys. Res. B 101 (1995) 499-510                                       
// Polarization of the outcoming photon:  Beam test of a prototype detector array for the PoGO astronomical hard X-ray/soft gamma-ray polarimeter
//                                        T. Mizuno et al., Nucl. Instr. Meth. in Phys. Res. A 540 (2005) 158-168                                        

  if (verboseLevel > 3)
    G4cout << "Calling G4LivermorePolarizedRayleighModel::Initialise()" << G4endl;

  if (crossSectionHandler)
  {
    crossSectionHandler->Clear();
    delete crossSectionHandler;
  }
  
  // Read data files for all materials

  crossSectionHandler = new G4CrossSectionHandler;
  crossSectionHandler->Clear();
  G4String crossSectionFile = "rayl/re-cs-";
  crossSectionHandler->LoadData(crossSectionFile);

  G4VDataSetAlgorithm* ffInterpolation = new G4LogLogInterpolation;
  G4String formFactorFile = "rayl/re-ff-";
  formFactorData = new G4CompositeEMDataSet(ffInterpolation,1.,1.);
  formFactorData->LoadData(formFactorFile);

  InitialiseElementSelectors(particle,cuts);

  //
  if (verboseLevel > 2) 
    G4cout << "Loaded cross section files for Livermore Polarized Rayleigh model" << G4endl;

  InitialiseElementSelectors(particle,cuts);

  if (verboseLevel > 0) { 
    G4cout << "Livermore Polarized Rayleigh model is initialized " << G4endl
         << "Energy range: "
         << LowEnergyLimit() / eV << " eV - "
         << HighEnergyLimit() / GeV << " GeV"
         << G4endl;
	 }

  if(isInitialised) return;
  fParticleChange = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedRayleighModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double GammaEnergy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
  if (verboseLevel > 3)
    G4cout << "Calling CrossSectionPerAtom() of G4LivermorePolarizedRayleighModel" << G4endl;

  if (GammaEnergy < lowEnergyLimit || GammaEnergy > highEnergyLimit) return 0.0;

  G4double cs = crossSectionHandler->FindValue(G4int(Z), GammaEnergy);
  return cs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedRayleighModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
					      const G4MaterialCutsCouple* couple,
					      const G4DynamicParticle* aDynamicGamma,
					      G4double,
					      G4double)
{
  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4LivermorePolarizedRayleighModel" << G4endl;

  G4double photonEnergy0 = aDynamicGamma->GetKineticEnergy();
  
  if (photonEnergy0 <= lowEnergyLimit)
  {
      fParticleChange->ProposeTrackStatus(fStopAndKill);
      fParticleChange->SetProposedKineticEnergy(0.);
      fParticleChange->ProposeLocalEnergyDeposit(photonEnergy0);
      return ;
  }

  G4ParticleMomentum photonDirection0 = aDynamicGamma->GetMomentumDirection();

  // Select randomly one element in the current material
  // G4int Z = crossSectionHandler->SelectRandomAtom(couple,photonEnergy0);
  const G4ParticleDefinition* particle =  aDynamicGamma->GetDefinition();
  const G4Element* elm = SelectRandomAtom(couple,particle,photonEnergy0);
  G4int Z = (G4int)elm->GetZ();

  G4double outcomingPhotonCosTheta = GenerateCosTheta(photonEnergy0, Z);
  G4double outcomingPhotonPhi = GeneratePhi(outcomingPhotonCosTheta);
  G4double beta=GeneratePolarizationAngle();
 
  // incomingPhoton reference frame:
  // z = versor parallel to the incomingPhotonDirection
  // x = versor parallel to the incomingPhotonPolarization
  // y = defined as z^x
 
  // outgoingPhoton reference frame:
  // z' = versor parallel to the outgoingPhotonDirection
  // x' = defined as x-x*z'z' normalized
  // y' = defined as z'^x'
 
  G4ThreeVector z(aDynamicGamma->GetMomentumDirection().unit()); 
  G4ThreeVector x(GetPhotonPolarization(*aDynamicGamma));
  G4ThreeVector y(z.cross(x));
 
  // z' = std::cos(phi)*std::sin(theta) x + std::sin(phi)*std::sin(theta) y + std::cos(theta) z
  G4double xDir;
  G4double yDir;
  G4double zDir;
  zDir=outcomingPhotonCosTheta;
  xDir=std::sqrt(1-outcomingPhotonCosTheta*outcomingPhotonCosTheta);
  yDir=xDir;
  xDir*=std::cos(outcomingPhotonPhi);
  yDir*=std::sin(outcomingPhotonPhi);
 
  G4ThreeVector zPrime((xDir*x + yDir*y + zDir*z).unit());
  G4ThreeVector xPrime(x.perpPart(zPrime).unit());
  G4ThreeVector yPrime(zPrime.cross(xPrime));
 
  // outgoingPhotonPolarization is directed as x' std::cos(beta) + y' std::sin(beta)
  G4ThreeVector outcomingPhotonPolarization(xPrime*std::cos(beta) + yPrime*std::sin(beta));
 
  fParticleChange->ProposeMomentumDirection(zPrime);
  fParticleChange->ProposePolarization(outcomingPhotonPolarization);
  fParticleChange->SetProposedKineticEnergy(photonEnergy0); 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedRayleighModel::GenerateCosTheta(G4double incomingPhotonEnergy, G4int zAtom) const
{
  //  d sigma                                                                    k0
  // --------- =  r0^2 * pi * F^2(x, Z) * ( 2 - sin^2 theta) * std::sin (theta), x = ---- std::sin(theta/2)
  //  d theta                                                                    hc
 
  //  d sigma                                             k0          1 - y
  // --------- = r0^2 * pi * F^2(x, Z) * ( 1 + y^2), x = ---- std::sqrt ( ------- ), y = std::cos(theta)
  //    d y                                               hc            2

  //              Z
  // F(x, Z) ~ --------
  //            a + bx
  //
  // The time to exit from the outer loop grows as ~ k0
  // On pcgeant2 the time is ~ 1 s for k0 ~ 1 MeV on the oxygen element. A 100 GeV
  // event will take ~ 10 hours.
  //
  // On the avarage the inner loop does 1.5 iterations before exiting
 
  const G4double xFactor = (incomingPhotonEnergy*cm)/(h_Planck*c_light);
  //const G4VEMDataSet * formFactorData = GetScatterFunctionData();

  G4double cosTheta;
  G4double fCosTheta;
  G4double x;
  G4double fValue;

  if (incomingPhotonEnergy > 5.*MeV)
  {
    cosTheta = 1.;
  }
  else
  {
    do
    {
      do
	{
	  cosTheta = 2.*G4UniformRand()-1.;
	  fCosTheta = (1.+cosTheta*cosTheta)/2.;
	}
      while (fCosTheta < G4UniformRand());
  
      x = xFactor*std::sqrt((1.-cosTheta)/2.);
  
      if (x > 1.e+005)
	fValue = formFactorData->FindValue(x, zAtom-1);
      else
	fValue = formFactorData->FindValue(0., zAtom-1);
   
      fValue/=zAtom;
      fValue*=fValue;
    }
    while(fValue < G4UniformRand());
  }

  return cosTheta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedRayleighModel::GeneratePhi(G4double cosTheta) const
{
  //  d sigma
  // --------- = alpha * ( 1 - sin^2 (theta) * cos^2 (phi) )
  //   d phi
 
  // On the average the loop takes no more than 2 iterations before exiting 

  G4double phi;
  G4double cosPhi;
  G4double phiProbability;
  G4double sin2Theta;
 
  sin2Theta=1.-cosTheta*cosTheta;
 
  do
    {
      phi = twopi * G4UniformRand();
      cosPhi = std::cos(phi);
      phiProbability= 1. - sin2Theta*cosPhi*cosPhi;
    }
  while (phiProbability < G4UniformRand());
 
  return phi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedRayleighModel::GeneratePolarizationAngle(void) const
{
  // Rayleigh polarization is always on the x' direction

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4LivermorePolarizedRayleighModel::GetPhotonPolarization(const G4DynamicParticle&  photon)
{

// SI - From G4VLowEnergyDiscretePhotonProcess.cc
 
  G4ThreeVector photonMomentumDirection;
  G4ThreeVector photonPolarization;

  photonPolarization = photon.GetPolarization(); 
  photonMomentumDirection = photon.GetMomentumDirection();

  if ((!photonPolarization.isOrthogonal(photonMomentumDirection, 1e-6)) || photonPolarization.mag()==0.)
    {
      // if |photonPolarization|==0. or |photonPolarization * photonDirection0| > 1e-6 * |photonPolarization ^ photonDirection0|
      // then polarization is choosen randomly.
  
      G4ThreeVector e1(photonMomentumDirection.orthogonal().unit());
      G4ThreeVector e2(photonMomentumDirection.cross(e1).unit());
  
      G4double angle(G4UniformRand() * twopi);
  
      e1*=std::cos(angle);
      e2*=std::sin(angle);
  
      photonPolarization=e1+e2;
    }
  else if (photonPolarization.howOrthogonal(photonMomentumDirection) != 0.)
    {
      // if |photonPolarization * photonDirection0| != 0.
      // then polarization is made orthonormal;
  
      photonPolarization=photonPolarization.perpPart(photonMomentumDirection);
    }
 
  return photonPolarization.unit();
}

