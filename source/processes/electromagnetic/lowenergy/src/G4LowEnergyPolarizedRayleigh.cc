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
// $Id: G4LowEnergyPolarizedRayleigh.cc,v 1.10 2009-06-11 15:47:08 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------
//
// File name:     G4LowEnergyPolarizedRayleigh.cc
//
// Author:        Capra Riccardo
//
// Creation date: May 2005
//
// History:
// -----------
// 02 May 2005  R. Capra         1st implementation
// 20 May 2005  MGP              Changed name of a local variable hiding 
//                               a data member of the base class
//
//----------------------------------------------------------------
      

#include "G4LowEnergyPolarizedRayleigh.hh"

#include "G4LogLogInterpolation.hh"
#include "G4VCrossSectionHandler.hh"
#include "G4VEMDataSet.hh"
#include "globals.hh" //G4Exception


G4LowEnergyPolarizedRayleigh::G4LowEnergyPolarizedRayleigh(const G4String& processName)
  :
  G4VLowEnergyDiscretePhotonProcess(processName, "rayl/re-cs-", "rayl/re-ff-", new G4LogLogInterpolation, 250*eV, 100*GeV),
  intrinsicLowEnergyLimit(10*eV),
  intrinsicHighEnergyLimit(100*GeV)
{
  if (GetLowEnergyLimit() < intrinsicLowEnergyLimit || 
      GetHighEnergyLimit() > intrinsicHighEnergyLimit)
    G4Exception("G4LowEnergyPolarizedRayleigh::G4LowEnergyPolarizedRayleigh - Energy limit outside intrinsic process validity range");
 
 
   G4cout << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "   The class G4LowEnergyPolarizedRayleigh is NOT SUPPORTED ANYMORE. " << G4endl;
   G4cout << "   It will be REMOVED with the next major release of Geant4. " << G4endl;
   G4cout << "   Please consult: https://twiki.cern.ch/twiki/bin/view/Geant4/LoweProcesses" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << G4endl;
}


G4VParticleChange* G4LowEnergyPolarizedRayleigh::PostStepDoIt(const G4Track&  aTrack, const G4Step&  aStep)
{
  // aParticleChange comes from G4VProcess
  aParticleChange.Initialize(aTrack);
 
  const G4DynamicParticle* incomingPhoton = aTrack.GetDynamicParticle();
  G4double incomingPhotonEnergy = incomingPhoton->GetKineticEnergy();
 
  if (incomingPhotonEnergy <= GetLowEnergyLimit())
    {
      aParticleChange.ProposeTrackStatus(fStopAndKill);
      aParticleChange.ProposeEnergy(0.);
      aParticleChange.ProposeLocalEnergyDeposit(incomingPhotonEnergy);
  
      return G4VLowEnergyDiscretePhotonProcess::PostStepDoIt(aTrack, aStep);
    }

  const G4VCrossSectionHandler* crossSectionHandle = GetCrossSectionHandler();
  const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
  G4int zAtom = crossSectionHandle->SelectRandomAtom(couple, incomingPhotonEnergy);

  G4double outcomingPhotonCosTheta = GenerateCosTheta(incomingPhotonEnergy, zAtom);
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
 
  G4ThreeVector z(incomingPhoton->GetMomentumDirection().unit()); 
  G4ThreeVector x(GetPhotonPolarization(*incomingPhoton));
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
 
  aParticleChange.ProposeEnergy(incomingPhotonEnergy);
  aParticleChange.ProposeMomentumDirection(zPrime);
  aParticleChange.ProposePolarization(outcomingPhotonPolarization);
  aParticleChange.SetNumberOfSecondaries(0);

  // returns aParticleChange though pParticleChange and G4VProcess::PostStepDoIt
  return G4VLowEnergyDiscretePhotonProcess::PostStepDoIt(aTrack, aStep);
}




G4double G4LowEnergyPolarizedRayleigh::GenerateCosTheta(G4double incomingPhotonEnergy, G4int zAtom) const
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
  const G4VEMDataSet * formFactorData = GetScatterFunctionData();

  G4double cosTheta;
  G4double fCosTheta;
  G4double x;
  G4double fValue;

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

  return cosTheta;
}



G4double G4LowEnergyPolarizedRayleigh::GeneratePhi(G4double cosTheta) const
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





G4double G4LowEnergyPolarizedRayleigh::GeneratePolarizationAngle(void) const
{
  // Rayleigh polarization is always on the x' direction

  return 0;
}
