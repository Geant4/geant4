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
/// \file exoticphysics/phonon/src/XPhononDownconversionProcess.cc
/// \brief Implementation of the XPhononDownconversionProcess class
//
// $Id$
//

#include "XPhononDownconversionProcess.hh"

#include "G4Step.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"

#include "XTPhononFast.hh"
#include "XTPhononSlow.hh"
#include "XLPhonon.hh"
#include <cmath>

#include "XPhononTrackInformation.hh"

#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4SystemOfUnits.hh"
#include "XLatticeManager3.hh"

#include "G4PhysicalConstants.hh"



XPhononDownconversionProcess::XPhononDownconversionProcess(const G4String& aName)
: G4VDiscreteProcess(aName)
{
   if (verboseLevel>1) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
   Lattice=0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XPhononDownconversionProcess::~XPhononDownconversionProcess()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XPhononDownconversionProcess::XPhononDownconversionProcess(XPhononDownconversionProcess& right)
: G4VDiscreteProcess(right)
{;}
 
G4double 
  XPhononDownconversionProcess::GetMeanFreePath( 
       const G4Track& aTrack, G4double /*previousStepSize*/, G4ForceCondition* condition  )
{
  //Determines mean free path for longitudinal phonons to split

  //Get pointer to lattice manager singleton and use it to find the
  //XPhysicalLattice object for current volume
  XLatticeManager3* LM = XLatticeManager3::GetXLatticeManager();
  Lattice = LM->GetXPhysicalLattice(aTrack.GetVolume());

  G4double A=Lattice->GetAnhDecConstant();
  G4double h=6.626068e-34*m2*kg/s; //Schroedinger's constant
  G4double E= aTrack.GetKineticEnergy();
  
  //Calculate mean free path for anh. decay
  G4double mfp = 1/((E/h)*(E/h)*(E/h)*(E/h)*(E/h)*A)*aTrack.GetVelocity();
  
  *condition = NotForced;
  return mfp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




G4VParticleChange*
  XPhononDownconversionProcess::PostStepDoIt( const G4Track& aTrack,
                                 const G4Step&)
{

  aParticleChange.Initialize(aTrack);

  //Get pointer to lattice manager singleton and use it to find the
  //XPhysicalLattice object for current volume
  XLatticeManager3* LM = XLatticeManager3::GetXLatticeManager();
  Lattice = LM->GetXPhysicalLattice(aTrack.GetVolume());

  //Obtain dynamical constants from this volume's lattice
  fBeta=Lattice->GetBeta();
  fGamma=Lattice->GetGamma();
  fLambda=Lattice->GetLambda();
  fMu=Lattice->GetMu();

  //Destroy the parent phonon and create the daughter phonons.
  //74% chance that daughter phonons are both transverse: call MakeTTSecondaries 
  //26% Transverse and Longitudinal: call MakeLTSecondaries
  if(G4UniformRand()>0.740) MakeLTSecondaries(aTrack); else MakeTTSecondaries(aTrack);
  aParticleChange.ProposeEnergy(0.);
  aParticleChange.ProposeTrackStatus(fStopAndKill);    
       
  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



G4bool XPhononDownconversionProcess::IsApplicable(const G4ParticleDefinition& aPD)
{
  //Only L-phonons decay
  return (&aPD==XLPhonon::PhononDefinition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


inline double XPhononDownconversionProcess::GetLTDecayProb(double d, double x){
  //probability density of energy distribution of L'-phonon in L->L'+T process
  //d=delta= ratio of group velocities vl/vt and x is the fraction of energy in the longitudinal mode, i.e. x=EL'/EL
  return (1/(x*x))*(1-x*x)*(1-x*x)*((1+x)*(1+x)-d*d*((1-x)*(1-x)))*(1+x*x-d*d*(1-x)*(1-x))*(1+x*x-d*d*(1-x)*(1-x));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


inline double XPhononDownconversionProcess::GetTTDecayProb(double d, double x){  
  //probability density of energy distribution of T-phonon in L->T+T process
  
  //dynamic constants from Tamura, PRL31, 1985
  G4double A = 0.5*(1-d*d)*(fBeta+fLambda+(1+d*d)*(fGamma+fMu));
  G4double B = fBeta+fLambda+2*d*d*(fGamma+fMu);
  G4double C = fBeta + fLambda + 2*(fGamma+fMu);
  G4double D = (1-d*d)*(2*fBeta+4*fGamma+fLambda+3*fMu);

  return (A+B*d*x-B*x*x)*(A+B*d*x-B*x*x)+(C*x*(d-x)-D/(d-x)*(x-d-(1-d*d)/(4*x)))*(C*x*(d-x)-D/(d-x)*(x-d-(1-d*d)/(4*x)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


inline double XPhononDownconversionProcess::MakeLDeviation(double d, double x){
  //change in L'-phonon propagation direction after decay

  return std::acos((1+(x*x)-((d*d)*(1-x)*(1-x)))/(2*x));
  //return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


inline double XPhononDownconversionProcess::MakeTDeviation(double d, double x){
  //change in T-phonon propagation direction after decay (L->L+T process)
  
  return std::acos((1-x*x+d*d*(1-x)*(1-x))/(2*d*(1-x)));
  //return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


inline double XPhononDownconversionProcess::MakeTTDeviation(double d, double x){
  //change in T-phonon propagation direction after decay (L->T+T process)

  return std::acos((1-d*d*(1-x)*(1-x)+d*d*x*x)/(2*d*x));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XPhononDownconversionProcess::MakeTTSecondaries(const G4Track& aTrack){
  //Generate daughter phonons from L->T+T process
   
  //d is the velocity ratio vL/vT
  G4double d=1.6338;
  G4double upperBound=(1+(1/d))/2;
  G4double lowerBound=(1-(1/d))/2;

  //Use MC method to generate point from distribution:
  //if a random point on the energy-probability plane is
  //smaller that the curve of the probability density,
  //then accept that point.
  //x=fraction of parent phonon energy in first T phonon
  G4double x = d*(G4UniformRand()*(upperBound-lowerBound)+lowerBound);
  G4double p = 1.5*G4UniformRand();
  while(!(p<GetTTDecayProb(d, x))){
    x=d*(G4UniformRand()*(upperBound-lowerBound)+lowerBound);
    p=1.5*G4UniformRand(); 
  }
  x=x/d;
  
  //using energy fraction x to calculate daughter phonon directions
  G4double theta1=MakeTTDeviation(d, (x));
  G4double theta2=MakeTTDeviation(d, (1-x));
  G4ThreeVector dir1=((XPhononTrackInformation*)(aTrack.GetUserInformation()))->GetK();
  G4ThreeVector dir2=dir1;
  G4ThreeVector ran = G4RandomDirection();
  
  //while(dir1.cross(ran)==0) ran=G4RandomDirection();
  //dir1 = dir1.rotate(dir1.cross(ran),theta1);
  //dir2 = dir2.rotate(dir2.cross(ran),-theta2);
  G4double ph=G4UniformRand()*2*pi;
  dir1 = dir1.rotate(dir1.orthogonal(),theta1).rotate(dir1, ph);
  dir2 = dir2.rotate(dir2.orthogonal(),-theta2).rotate(dir2,ph);

  aParticleChange.SetNumberOfSecondaries(2);
  G4double E=aTrack.GetKineticEnergy();
  G4Track* sec1;
  G4Track* sec2;
  
  G4double probST = Lattice->GetSTDOS()/(Lattice->GetSTDOS()+Lattice->GetFTDOS());

 //First secondary:Make FT or ST phonon, probability density is funciton of eqn of state
  int polarization1;
  if(G4UniformRand()<probST){ // DOS_slow / (DOS_slow+DOS_fast) = 0.59345 according to ModeDensity.m
    polarization1 = 1;
    sec1 = new G4Track(new G4DynamicParticle(XTPhononSlow::PhononDefinition(),Lattice->MapKtoVDir(1,dir1), x*E),aTrack.GetGlobalTime(), aTrack.GetPosition() );
  }else{
    polarization1 = 2;
    sec1 = new G4Track(new G4DynamicParticle(XTPhononFast::PhononDefinition(),Lattice->MapKtoVDir(2,dir1), x*E),aTrack.GetGlobalTime(), aTrack.GetPosition() );
  }

 //Second secondary:Make FT or ST phonon, probability density is funciton of eqn of state
    int polarization2;
  if(G4UniformRand()<probST){ // DOS_slow / (DOS_slow+DOS_fast) = 0.59345 according to ModeDensity.m
    polarization2 = 1;
    sec2 = new G4Track(new G4DynamicParticle(XTPhononSlow::PhononDefinition(),Lattice->MapKtoVDir(1,dir2), (1-x)*E),aTrack.GetGlobalTime(), aTrack.GetPosition() );
  }else{
    polarization2 = 2;
    sec2 = new G4Track(new G4DynamicParticle(XTPhononFast::PhononDefinition(),Lattice->MapKtoVDir(2,dir2), (1-x)*E),aTrack.GetGlobalTime(), aTrack.GetPosition() );
  }

  //Set the k-vectors for the two secondaries and add them to the process
  sec1->SetUserInformation(new XPhononTrackInformation(dir1));
  sec1->SetVelocity(Lattice->MapKtoV(polarization1, dir1)*m/s);
  sec1->UseGivenVelocity(true);

  sec2->SetUserInformation(new XPhononTrackInformation(dir2));
  sec2->SetVelocity(Lattice->MapKtoV(polarization2, dir2)*m/s);
  sec2->UseGivenVelocity(true);

  aParticleChange.AddSecondary(sec1);
  aParticleChange.AddSecondary(sec2);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XPhononDownconversionProcess::MakeLTSecondaries(const G4Track& aTrack){

  //d is the velocity ratio vL/v
  G4double d=1.6338;
  G4double upperBound=1;
  G4double lowerBound=(d-1)/(d+1);
  
  //Use MC method to generate point from distribution:
  //if a random point on the energy-probability plane is
  //smaller that the curve of the probability density,
  //then accept that point.
  //x=fraction of parent phonon energy in L phonon
  G4double x =(G4UniformRand()*(upperBound-lowerBound)+lowerBound);
  G4double p = 4.0*G4UniformRand();
  while(!(p<GetLTDecayProb(d, x))){
    x=(G4UniformRand()*(upperBound-lowerBound)+lowerBound);
    p=4.0*G4UniformRand(); //4.0 is about the max in the probability density function
  }

  //using energy fraction x to calculate daughter phonon directions
  G4double thetaL=MakeLDeviation(d, x);
  G4double thetaT=MakeTDeviation(d, x);;
  G4ThreeVector dir1=((XPhononTrackInformation*)(aTrack.GetUserInformation()))->GetK();
  G4ThreeVector dir2=dir1;

  G4double ph=G4UniformRand()*2*pi;
  dir1 = dir1.rotate(dir1.orthogonal(),thetaL).rotate(dir1, ph);
  dir2 = dir2.rotate(dir2.orthogonal(),-thetaT).rotate(dir2,ph);



  aParticleChange.SetNumberOfSecondaries(2);
  G4double E=aTrack.GetKineticEnergy();

  G4Track* sec1 = new G4Track(new G4DynamicParticle(XLPhonon::PhononDefinition(),Lattice->MapKtoVDir(0,dir1), x*E),aTrack.GetGlobalTime(), aTrack.GetPosition() );

  G4Track* sec2;

  //Make FT or ST phonon, probability density is funciton of eqn of state
  G4double probST = Lattice->GetSTDOS()/(Lattice->GetSTDOS()+Lattice->GetFTDOS());
  int polarization;
  if(G4UniformRand()<probST){ // DOS_slow / (DOS_slow+DOS_fast) = 0.59345 according to ModeDensity.m
    
    sec2 = new G4Track(new G4DynamicParticle(XTPhononSlow::PhononDefinition(),Lattice->MapKtoVDir(1,dir2), (1-x)*E),aTrack.GetGlobalTime(), aTrack.GetPosition() );
    polarization = 1;
  
  }else{
    
    sec2 = new G4Track(new G4DynamicParticle(XTPhononFast::PhononDefinition(),Lattice->MapKtoVDir(2,dir2), (1-x)*E),aTrack.GetGlobalTime(), aTrack.GetPosition() );
    polarization = 2;

  }

  sec1->SetUserInformation(new XPhononTrackInformation(dir1));
  sec1->SetVelocity(Lattice->MapKtoV(0, dir1)*m/s);
  sec1->UseGivenVelocity(true);
  
  sec2->SetUserInformation(new XPhononTrackInformation(dir2));
  sec2->SetVelocity(Lattice->MapKtoV(polarization, dir2)*m/s);
  sec2->UseGivenVelocity(true);

  aParticleChange.AddSecondary(sec1);
  aParticleChange.AddSecondary(sec2); 

    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

