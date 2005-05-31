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
//
// $Id: G4DNAElectronElasticScatteringInWater.cc,v 1.1 2005-05-31 09:58:40 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4DNAElectronElasticScatteringInWater.hh"

#include "G4Electron.hh"
 
                                         G4DNAElectronElasticScatteringInWater :: G4DNAElectronElasticScatteringInWater(const G4String & name)
:
 G4VDNAProcessInWater(name),
 lowEnergyLimit(7*eV),
 highEnergyLimit(10*keV),
 intrinsicLowEnergyLimit(7*eV),
 intrinsicHighEnergyLimit(10*keV)
{
 if (lowEnergyLimit < intrinsicLowEnergyLimit ||
     highEnergyLimit > intrinsicHighEnergyLimit)
  G4Exception("G4DNAElectronElasticScatteringInWater::G4DNAElectronElasticScatteringInWater: Energy outside intrinsic process validity range");
}

G4VParticleChange *                      G4DNAElectronElasticScatteringInWater :: PostStepDoIt(const G4Track & aTrack, const G4Step & aStep)
{
 aParticleChange.Initialize(aTrack);
 
 G4double k;
 k=aTrack.GetDynamicParticle()->GetKineticEnergy();

 if (k <= lowEnergyLimit)
 {
  aParticleChange.ProposeTrackStatus(fStopAndKill);
  aParticleChange.ProposeEnergy(0.);
  aParticleChange.ProposeLocalEnergyDeposit(k);
  
  return G4VDNAProcessInWater::PostStepDoIt(aTrack, aStep);
 }
 
 G4double cosTheta;
 const G4int z(10); // H2O number of electrons
 
 if (k > 200*eV)
  cosTheta=GenerateCosThetaElasticEmfietzoglou(k, z);
 else
  cosTheta=GenerateCosThetaElasticBrenner(k);

 G4double phi;
 phi=2*pi*G4UniformRand();
 
 G4ThreeVector zVers(aTrack.GetDynamicParticle()->GetMomentumDirection());
 G4ThreeVector xVers(zVers.orthogonal());
 G4ThreeVector yVers(zVers.cross(xVers));
 
 G4double xDir;
 G4double yDir;
 
 xDir=std::sqrt(1-cosTheta*cosTheta);
 yDir=xDir;
 xDir*=cos(phi);
 yDir*=sin(phi);
 
 G4ThreeVector zPrimeVers((xDir*xVers + yDir*yVers + cosTheta*zVers).unit());

 aParticleChange.ProposeEnergy(k);
 aParticleChange.ProposeMomentumDirection(zPrimeVers);
 aParticleChange.SetNumberOfSecondaries(0);
 
 return G4VDNAProcessInWater::PostStepDoIt(aTrack, aStep);
}

G4bool                                   G4DNAElectronElasticScatteringInWater :: IsApplicable(const G4ParticleDefinition & aParticleDefinition)
{
 return (&aParticleDefinition) == G4Electron::Electron(); 
}

G4double                                 G4DNAElectronElasticScatteringInWater :: GetMeanFreePath(const G4Track & aTrack, G4double /*previousStepSize*/, G4ForceCondition * /*condition*/)
{
 //                pi * sigma_Ruth(K)
 // sigma_el(K) = --------------------
 //                n(K) * (n(K) + 1)
 //
 // Where K is the electron non-relativistic kinetic energy
 // Cross section per water molecule
 
 const G4int z(10); // H2O number of electrons
 
 G4double k;
 k=aTrack.GetDynamicParticle()->GetKineticEnergy();
 
 if (k < lowEnergyLimit)
  k=lowEnergyLimit;
 
 if (k > highEnergyLimit)
  k=highEnergyLimit;
 
 G4double n;
 n=ScreeningFactor(k, z);
 
 G4double sigma_Ruth;
 sigma_Ruth=RutherfordTotalCrossSection(k, z);
 
 G4double sigma_el;
 sigma_el=(pi*sigma_Ruth)/(n*(n+1.));
 
 // We suppose we are in water, one of the elements must be oxygen
 G4Material * theMaterial(aTrack.GetMaterial());

 size_t i(theMaterial->GetNumberOfElements());
 
 while (i>0)
 {
  i--;
  
  const G4Element * element(theMaterial->GetElement(i));
  
  if (element->GetZ()==8.)
  {
   // Number of oxigens per volume = number of water molecules per volume
   G4double density;
   density=theMaterial->GetAtomicNumDensityVector()[i];
  
   return 1./(density*sigma_el);
  }
 }
 
 G4String message;
 message="G4DNAElectronElasticScatteringInWater::GetMeanFreePath - ";
 message+=theMaterial->GetName();
 message+=" is not a water material";
 
 G4Exception(message);
 
 return DBL_MAX;
}
  
G4double                                 G4DNAElectronElasticScatteringInWater :: RutherfordTotalCrossSection(G4double k, G4int z)
{
 //                                  e^4         /      K + m_e c^2      \^2
 // sigma_Ruth(K) = Z (Z+1) -------------------- | --------------------- |
 //                          (4 pi epsilon_0)^2  \  K * (K + 2 m_e c^2)  /
 //
 // Where K is the electron non-relativistic kinetic energy
 // 
 // Nucl. Instr. Meth. 155 (1978) 145-156
 
 G4double length;
 length=(e_squared*(k+electron_mass_c2))/(4*pi*epsilon0*k*(k+2*electron_mass_c2));
 
 return static_cast<G4double>(z*(z+1))*length*length;
}

G4double                                 G4DNAElectronElasticScatteringInWater :: ScreeningFactor(G4double k, G4int z)
{
 //
 //         alpha_1 + beta_1 ln(K/eV)   constK Z^(2/3)
 // n(T) = -------------------------- -----------------
 //              K/(m_e c^2)            2 + K/(m_e c^2)
 //
 // Where K is the electron non-relativistic kinetic energy
 //
 // n(T) > 0 for T < ~ 400 MeV
 // 
 // Nucl. Instr. Meth. 155 (1978) 145-156

 const G4double alpha_1(1.64);
 const G4double beta_1(-0.0825);
 const G4double constK(1.7E-5);
 
 G4double numerator;
 numerator=(alpha_1+beta_1*std::log(k/eV))*constK*std::pow(static_cast<double>(z), 2./3.);

 k/=electron_mass_c2;

 G4double denominator;
 denominator=k*(2+k);
 
 return numerator/denominator;
}

G4double                                 G4DNAElectronElasticScatteringInWater :: GenerateCosThetaElasticEmfietzoglou(G4double k, G4int z)
{
 //  d sigma_el                sigma_Ruth(K)
 // ------------ (K) ~ -----------------------------  
 //   d Omega           (1 + 2 n(K) - cos(theta))^2
 //
 // We extract cos(theta) distributed as (1 + 2 n(K) - cos(theta))^-2
 //
 // Maximum is for theta=0: 1/(4 n(K)^2) (When n(K) is positive, that is always satisfied within the validity of the process)
 //
 // Phys. Med. Biol. 45 (2000) 3171-3194
 
 G4double n;
 n=ScreeningFactor(k, z);

 G4double oneOverMax;
 oneOverMax=(4.*n*n);
 
 G4double cosTheta;
 G4double fCosTheta;
 
 do
 {
  cosTheta = 2.*G4UniformRand()-1.;
  fCosTheta = (1 + 2.*n - cosTheta);
  fCosTheta = oneOverMax/(fCosTheta*fCosTheta);
 }
 while (fCosTheta < G4UniformRand());
 
 return cosTheta;
}

G4double                                 G4DNAElectronElasticScatteringInWater :: GenerateCosThetaElasticBrenner(G4double k)
{
 //  d sigma_el                         1                                 beta(K)
 // ------------ (K) ~ --------------------------------- + ---------------------------------
 //   d Omega           (1 + 2 gamma(K) - cos(theta))^2     (1 + 2 delta(K) + cos(theta))^2
 //
 // Maximum is < 1/(4 gamma(K)^2) + beta(K)/(4 delta(K)^2)
 //
 // Phys. Med. Biol. 29 N.4 (1983) 443-447
 
 const G4double betaCoeff[5]         = { 7.51525,   -0.419122,    7.2017E-3, -4.646E-5,    1.02897E-7};  
 const G4double deltaCoeff[5]        = { 2.9612,    -0.26376,     4.307E-3,  -2.6895E-5,   5.83505E-8};
 const G4double gamma035_10Coeff[6]  = {-1.7013,    -1.48284,     0.6331,    -0.10911,     8.358E-3,  -2.388E-4};
 const G4double gamma10_100Coeff[5]  = {-3.32517,    0.10996,    -4.5255E-3,  5.8372E-5,  -2.4659E-7};
 const G4double gamma100_200Coeff[3] = { 2.4775E-2, -2.96264E-5, -1.20655E-7};

 // gamma(K), beta(K) and delta(K) are polynomials with coefficients for energy measured in eV
 k/=eV;
 
 G4double beta;
 beta=std::exp(CalulatePolynomial(k, betaCoeff, 5)); 

 G4double delta;
 delta=std::exp(CalulatePolynomial(k, deltaCoeff, 5)); 

 G4double gamma;
 
 if (k>100)
  gamma=CalulatePolynomial(k, gamma100_200Coeff, 3); // Only in this case it is not the exponent of the polynomial
 else if (k>10)
  gamma=std::exp(CalulatePolynomial(k, gamma10_100Coeff, 5));
 else
  gamma=std::exp(CalulatePolynomial(k, gamma035_10Coeff, 6));
  
 G4double oneOverMax;
 oneOverMax=1./(1./(4.*gamma*gamma) + beta/(4.*delta*delta));
 
 G4double cosTheta;
 G4double leftDenominator;
 G4double rightDenominator;
 G4double fCosTheta;
 
 do
 {
  cosTheta = 2.*G4UniformRand()-1.;
  leftDenominator = (1 + 2.*gamma - cosTheta);
  rightDenominator = (1 + 2*delta + cosTheta);
  fCosTheta = oneOverMax*(1./(leftDenominator*leftDenominator) + beta/(rightDenominator*rightDenominator));
 }
 while (fCosTheta < G4UniformRand());
 
 return cosTheta; 
}

G4double                                 G4DNAElectronElasticScatteringInWater :: CalulatePolynomial(G4double k, const G4double *vector, G4int size)
{
 // Sum_{i=0}^{size-1} vector_i k^i
 //
 // Phys. Med. Biol. 29 N.4 (1983) 443-447

 G4double result(0);

 while (size>0)
 {
  size--;
  
  result*=k;
  result+=vector[size];
 }
 
 return result;
}
