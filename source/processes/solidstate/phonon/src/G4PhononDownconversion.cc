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
/// \file processes/phonon/src/G4PhononDownconversion.cc
/// \brief Implementation of the G4PhononDownconversion class
//
// $Id: G4PhononDownconversion.cc 79766 2014-03-13 16:11:36Z gcosmo $
//
// 20131111  Add verbose output for MFP calculation
// 20131115  Initialize data buffers in ctor

#include "G4PhononDownconversion.hh"
#include "G4LatticePhysical.hh"
#include "G4PhononLong.hh"
#include "G4PhononPolarization.hh"
#include "G4PhononTrackMap.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include <cmath>



G4PhononDownconversion::G4PhononDownconversion(const G4String& aName)
  : G4VPhononProcess(aName), fBeta(0.), fGamma(0.), fLambda(0.), fMu(0.) {;}

G4PhononDownconversion::~G4PhononDownconversion() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PhononDownconversion::GetMeanFreePath(const G4Track& aTrack,
						 G4double /*previousStepSize*/,
						 G4ForceCondition* condition) {
  //Determines mean free path for longitudinal phonons to split
  G4double A = theLattice->GetAnhDecConstant();
  G4double Eoverh = aTrack.GetKineticEnergy()/h_Planck;
  
  //Calculate mean free path for anh. decay
  G4double mfp = aTrack.GetVelocity()/(Eoverh*Eoverh*Eoverh*Eoverh*Eoverh*A);

  if (verboseLevel > 1)
    G4cout << "G4PhononDownconversion::GetMeanFreePath = " << mfp << G4endl;
  
  *condition = NotForced;
  return mfp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4VParticleChange* G4PhononDownconversion::PostStepDoIt( const G4Track& aTrack,
							 const G4Step&) {
  aParticleChange.Initialize(aTrack);

  //Obtain dynamical constants from this volume's lattice
  fBeta=theLattice->GetBeta();
  fGamma=theLattice->GetGamma();
  fLambda=theLattice->GetLambda();
  fMu=theLattice->GetMu();

  //Destroy the parent phonon and create the daughter phonons.
  //74% chance that daughter phonons are both transverse
  //26% Transverse and Longitudinal
  if (G4UniformRand()>0.740) MakeLTSecondaries(aTrack);
  else MakeTTSecondaries(aTrack);

  aParticleChange.ProposeEnergy(0.);
  aParticleChange.ProposeTrackStatus(fStopAndKill);    
       
  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4PhononDownconversion::IsApplicable(const G4ParticleDefinition& aPD) {
  //Only L-phonons decay
  return (&aPD==G4PhononLong::PhononDefinition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//probability density of energy distribution of L'-phonon in L->L'+T process

G4double G4PhononDownconversion::GetLTDecayProb(G4double d, G4double x) const {
  //d=delta= ratio of group velocities vl/vt and x is the fraction of energy in the longitudinal mode, i.e. x=EL'/EL
  return (1/(x*x))*(1-x*x)*(1-x*x)*((1+x)*(1+x)-d*d*((1-x)*(1-x)))*(1+x*x-d*d*(1-x)*(1-x))*(1+x*x-d*d*(1-x)*(1-x));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//probability density of energy distribution of T-phonon in L->T+T process
G4double G4PhononDownconversion::GetTTDecayProb(G4double d, G4double x) const {  
  //dynamic constants from Tamura, PRL31, 1985
  G4double A = 0.5*(1-d*d)*(fBeta+fLambda+(1+d*d)*(fGamma+fMu));
  G4double B = fBeta+fLambda+2*d*d*(fGamma+fMu);
  G4double C = fBeta + fLambda + 2*(fGamma+fMu);
  G4double D = (1-d*d)*(2*fBeta+4*fGamma+fLambda+3*fMu);

  return (A+B*d*x-B*x*x)*(A+B*d*x-B*x*x)+(C*x*(d-x)-D/(d-x)*(x-d-(1-d*d)/(4*x)))*(C*x*(d-x)-D/(d-x)*(x-d-(1-d*d)/(4*x)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4double G4PhononDownconversion::MakeLDeviation(G4double d, G4double x) const {
  //change in L'-phonon propagation direction after decay

  return std::acos((1+(x*x)-((d*d)*(1-x)*(1-x)))/(2*x));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4double G4PhononDownconversion::MakeTDeviation(G4double d, G4double x) const {
  //change in T-phonon propagation direction after decay (L->L+T process)
  
  return std::acos((1-x*x+d*d*(1-x)*(1-x))/(2*d*(1-x)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4double G4PhononDownconversion::MakeTTDeviation(G4double d, G4double x) const {
  //change in T-phonon propagation direction after decay (L->T+T process)

  return std::acos((1-d*d*(1-x)*(1-x)+d*d*x*x)/(2*d*x));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


//Generate daughter phonons from L->T+T process
   
void G4PhononDownconversion::MakeTTSecondaries(const G4Track& aTrack) {
  //d is the velocity ratio vL/vT
  G4double d=1.6338;
  G4double upperBound=(1+(1/d))/2;
  G4double lowerBound=(1-(1/d))/2;

  //Use MC method to generate point from distribution:
  //if a random point on the energy-probability plane is
  //smaller that the curve of the probability density,
  //then accept that point.
  //x=fraction of parent phonon energy in first T phonon
  G4double x = G4UniformRand()*(upperBound-lowerBound) + lowerBound;
  G4double p = 1.5*G4UniformRand();
  while(p >= GetTTDecayProb(d, x*d)) {
    x = G4UniformRand()*(upperBound-lowerBound) + lowerBound;
    p = 1.5*G4UniformRand(); 
  }
  
  //using energy fraction x to calculate daughter phonon directions
  G4double theta1=MakeTTDeviation(d, x);
  G4double theta2=MakeTTDeviation(d, 1-x);
  G4ThreeVector dir1=trackKmap->GetK(aTrack);
  G4ThreeVector dir2=dir1;

  // FIXME:  These extra randoms change timing and causting outputs of example!
  G4ThreeVector ran = G4RandomDirection();	// FIXME: Drop this line
  
  G4double ph=G4UniformRand()*twopi;
  dir1 = dir1.rotate(dir1.orthogonal(),theta1).rotate(dir1, ph);
  dir2 = dir2.rotate(dir2.orthogonal(),-theta2).rotate(dir2,ph);

  G4double E=aTrack.GetKineticEnergy();
  G4double Esec1 = x*E, Esec2 = E-Esec1;

  // Make FT or ST phonon (0. means no longitudinal)
  G4int polarization1 = ChoosePolarization(0., theLattice->GetSTDOS(),
					   theLattice->GetFTDOS());

  // Make FT or ST phonon (0. means no longitudinal)
  G4int polarization2 = ChoosePolarization(0., theLattice->GetSTDOS(),
					   theLattice->GetFTDOS());

  // Construct the secondaries and set their wavevectors
  G4Track* sec1 = CreateSecondary(polarization1, dir1, Esec1);
  G4Track* sec2 = CreateSecondary(polarization2, dir2, Esec2);

  aParticleChange.SetNumberOfSecondaries(2);
  aParticleChange.AddSecondary(sec1);
  aParticleChange.AddSecondary(sec2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


//Generate daughter phonons from L->L'+T process
   
void G4PhononDownconversion::MakeLTSecondaries(const G4Track& aTrack) {
  //d is the velocity ratio vL/v
  G4double d=1.6338;
  G4double upperBound=1;
  G4double lowerBound=(d-1)/(d+1);
  
  //Use MC method to generate point from distribution:
  //if a random point on the energy-probability plane is
  //smaller that the curve of the probability density,
  //then accept that point.
  //x=fraction of parent phonon energy in L phonon
  G4double x = G4UniformRand()*(upperBound-lowerBound) + lowerBound;
  G4double p = 4.0*G4UniformRand();
  while(p >= GetLTDecayProb(d, x)) {
    x = G4UniformRand()*(upperBound-lowerBound) + lowerBound;
    p = 4.0*G4UniformRand(); 		     //4.0 is about the max in the PDF
  }

  //using energy fraction x to calculate daughter phonon directions
  G4double thetaL=MakeLDeviation(d, x);
  G4double thetaT=MakeTDeviation(d, x);		// FIXME:  Should be 1-x?
  G4ThreeVector dir1=trackKmap->GetK(aTrack);
  G4ThreeVector dir2=dir1;

  G4double ph=G4UniformRand()*twopi;
  dir1 = dir1.rotate(dir1.orthogonal(),thetaL).rotate(dir1, ph);
  dir2 = dir2.rotate(dir2.orthogonal(),-thetaT).rotate(dir2,ph);

  G4double E=aTrack.GetKineticEnergy();
  G4double Esec1 = x*E, Esec2 = E-Esec1;

  // First secondary is longitudnal
  G4int polarization1 = G4PhononPolarization::Long;

  // Make FT or ST phonon (0. means no longitudinal)
  G4int polarization2 = ChoosePolarization(0., theLattice->GetSTDOS(),
					   theLattice->GetFTDOS());

  // Construct the secondaries and set their wavevectors
  G4Track* sec1 = CreateSecondary(polarization1, dir1, Esec1);
  G4Track* sec2 = CreateSecondary(polarization2, dir2, Esec2);

  aParticleChange.SetNumberOfSecondaries(2);
  aParticleChange.AddSecondary(sec1);
  aParticleChange.AddSecondary(sec2); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

