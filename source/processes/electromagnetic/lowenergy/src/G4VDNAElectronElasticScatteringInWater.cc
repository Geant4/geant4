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
// $Id: G4VDNAElectronElasticScatteringInWater.cc,v 1.1 2005-06-02 15:02:54 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4VDNAElectronElasticScatteringInWater.hh"

#include "G4Electron.hh"
 
                                         G4VDNAElectronElasticScatteringInWater :: G4VDNAElectronElasticScatteringInWater(const G4String & name)
:
 G4VDNAProcessInWater(name)
{
}

G4VParticleChange *                      G4VDNAElectronElasticScatteringInWater :: PostStepDoIt(const G4Track & aTrack, const G4Step & aStep)
{
 aParticleChange.Initialize(aTrack);
 
 G4double k;
 k=aTrack.GetDynamicParticle()->GetKineticEnergy();

 G4double cosTheta;
 const G4int z(10); // H2O number of electrons
 
 cosTheta=RandomizeCosTheta(k, z);

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

G4bool                                   G4VDNAElectronElasticScatteringInWater :: IsApplicable(const G4ParticleDefinition & aParticleDefinition)
{
 return (&aParticleDefinition) == G4Electron::Electron(); 
}

G4double                                 G4VDNAElectronElasticScatteringInWater :: RutherfordTotalCrossSection(G4double k, G4int z) const
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

G4double                                 G4VDNAElectronElasticScatteringInWater :: ScreeningFactor(G4double k, G4int z) const
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

