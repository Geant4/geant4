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
//
// $Id: G4CrossSectionExcitationMillerGreen.cc,v 1.1 2007-05-02 17:20:36 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// $Id: G4CrossSectionExcitationMillerGreen.cc,v 1.1 2007-05-02 17:20:36 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Contact Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// Reference: TNS Geant4-DNA paper
// Reference for implementation model: NIM. 155, pp. 145-156, 1978

// History:
// -----------
// Date         Name              Modification
// 28 Apr 2007  M.G. Pia          Created in compliance with design described in TNS paper
//
// -------------------------------------------------------------------

// Class description:
// Geant4-DNA Cross total cross section for electron elastic scattering in water
// Further documentation available from http://www.ge.infn.it/geant4/dna

// -------------------------------------------------------------------


#include "G4CrossSectionExcitationMillerGreen.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"

#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Electron.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"

#include "Randomize.hh"

G4CrossSectionExcitationMillerGreen::G4CrossSectionExcitationMillerGreen()
{

  name = "CrossSectionExcitationMillerGreen";
  lowEnergyLimit = 7. * eV;
  highEnergyLimit = 10 * keV;

  // Energy correction factor for protons (to be added for alpha and other particles)
  kineticEnergyCorrection = 1.;

  // The following are for protons; to be generalized for alpha and other particles
  slaterEffectiveCharge[0] = 0.;
  slaterEffectiveCharge[1] = 0.;
  slaterEffectiveCharge[2] = 0.;
  sCoefficient[0] = 0.;
  sCoefficient[1] = 0.;
  sCoefficient[2] = 0.;

//  if (verboseLevel > 0)
//  {
//    G4cout << name << " is created " << G4endl
//     << "Energy range: "
//     << lowEnergyLimit / keV << " keV - "
//     << highEnergyLimit / GeV << " GeV"
//     << G4endl;
//  }
}


G4CrossSectionExcitationMillerGreen::~G4CrossSectionExcitationMillerGreen()
{ }
 

G4double G4CrossSectionExcitationMillerGreen::CrossSection(const G4Track& track)
{
  const G4DynamicParticle* particle = track.GetDynamicParticle();
  G4double k = particle->GetKineticEnergy();
  const G4ParticleDefinition* particleDef = track.GetDefinition();

  // G4Material* material = track.GetMaterial();

  // Assume that the material is water; proper algorithm to calculate correctly for any material to be inserted here
  G4int z = 10;

  // Take into account 5 excitation levels (D. Emfietzoglou et al., NIM B 193, pp. 71-78, 2002.
  G4int i = 5;

  G4double crossSection(0.);

  while (i>0)
    {
      i--;
      crossSection += PartialCrossSection(k,z,i,particleDef);
    }
  return crossSection;
 }


G4double G4CrossSectionExcitationMillerGreen::PartialCrossSection(G4double k, 
								  G4int z, 
								  G4int excitationLevel, 
								  const G4ParticleDefinition* particle)
{
   //                               ( ( z * aj ) ^ omegaj ) * ( t - ej ) ^ nu
  // sigma(t) = zEff^2 * sigma0 * --------------------------------------------
  //                               jj ^ ( omegaj + nu ) + t ^ ( omegaj + nu )
  //
  // where t is the kinetic energy corrected by Helium mass over proton mass for Helium ions
  //
  // zEff is:
  //  1 for protons
  //  2 for alpha++
  //  and  2 - c1 S_1s - c2 S_2s - c3 S_2p for alpha+ and He
  //
  // Dingfelder et al., RPC 59 p. 266 (2000) from Miller and Green (1973)
  
  const G4double sigma0(1.E+8 * barn);
  const G4double nu(1.);
  const G4double aj[]={876.*eV, 2084.* eV, 1373.*eV, 692.*eV, 900.*eV};
  const G4double jj[]={19820.*eV, 23490.*eV, 27770.*eV, 30830.*eV, 33080.*eV};
  const G4double omegaj[]={0.85, 0.88, 0.88, 0.78, 0.78};
  
  G4double tCorrected;
  tCorrected = k * kineticEnergyCorrection;

  G4double numerator;
  numerator = std::pow(z * aj[excitationLevel], omegaj[excitationLevel]) * 
    std::pow(tCorrected - EnergyConstant(excitationLevel), nu);

  G4double power;
  power = omegaj[excitationLevel] + nu;

  G4double denominator;
  denominator = std::pow(jj[excitationLevel], power) + std::pow(tCorrected, power);

  G4double zEff = particle->GetPDGCharge() / eplus + particle->GetLeptonNumber();
  
  zEff -= ( sCoefficient[0] * S_1s(k, EnergyConstant(excitationLevel), slaterEffectiveCharge[0], 1.) +
	    sCoefficient[1] * S_2s(k, EnergyConstant(excitationLevel), slaterEffectiveCharge[1], 2.) +
	    sCoefficient[2] * S_2p(k, EnergyConstant(excitationLevel), slaterEffectiveCharge[2], 2.) );

  G4double cross = sigma0 * zEff * zEff * numerator / denominator;

  return cross;
}


G4double G4CrossSectionExcitationMillerGreen::S_1s(G4double t, 
						   G4double energyTransferred, 
						   G4double slaterEffectiveCharge, 
						   G4double shellNumber)
{
  // 1 - e^(-2r) * ( 1 + 2 r + 2 r^2)
 
  G4double r = R(t, energyTransferred, slaterEffectiveCharge, shellNumber);
  G4double value = 1. - std::exp(-2 * r) * ( ( 2. * r + 2. ) * r + 1. );
  
  return value;
}


G4double G4CrossSectionExcitationMillerGreen::S_2s(G4double t, 
						   G4double energyTransferred, 
						   G4double slaterEffectiveCharge, 
						   G4double shellNumber)
{
  // 1 - e^(-2 r) * ( 1 + 2 r + 2 r^2 + 2 r^4)

  G4double r = R(t, energyTransferred, slaterEffectiveCharge, shellNumber);
  G4double value =  1. - std::exp(-2 * r) * (((2. * r * r + 2.) * r + 2.) * r + 1.);

  return value;
 
}

G4double G4CrossSectionExcitationMillerGreen::S_2p(G4double t, 
						   G4double energyTransferred, 
						   G4double slaterEffectiveCharge, 
						   G4double shellNumber)
{

  // 1 - e^(-2 r) * ( 1 + 2 r + 2 r^2 + 4/3 r^3 + 2/3 r^4)

  G4double r = R(t, energyTransferred, slaterEffectiveCharge, shellNumber);
  G4double value =  1. - std::exp(-2 * r) * (((( 2./3. * r + 4./3.) * r + 2.) * r + 2.) * r  + 1.);

  return value;
}


G4double G4CrossSectionExcitationMillerGreen::R(G4double t, 
						G4double energyTransferred, 
						G4double slaterEffectiveCharge, 
						G4double shellNumber) 
{
  // tElectron = m_electron / m_alpha * t

  // Hardcoded in Riccardo's implementation; to be corrected
  G4double tElectron = 0.511/3728. * t;
  G4double value = 2. * tElectron * slaterEffectiveCharge / (energyTransferred * shellNumber);
  
  return value;
}



 //G4int G4CrossSectionExcitationMillerGreen::RandomizePartialCrossSection(G4double k, G4int z)
 //{
  // Assume 5 excitation levels 

  //  G4int i = 5;
  // G4double value = 0.;
  //G4double values[5];
  
   //while (i>0)
     //{
   //i--;
    //   values[i] = PartialCrossSection(k,z,i,definition);
    //   value += values[i];
    // }
  
 //  value *= G4UniformRand();
  
   // i = 5;
  // while (i>0)
   //  {
     //  i--;
      
     //  if (values[i] > value) return i;
      
     //  value -= values[i];
  //   }
  // One should never end up here; next statement added to avoid compilation warning
  // Probably one should throw an exception if one ends up here
 //  return 0;
 //}

G4double G4CrossSectionExcitationMillerGreen::EnergyConstant(G4int excitationLevel)
{
  const G4double ej[]={8.17*eV, 10.13*eV, 11.31*eV, 12.91*eV, 14.50*eV};
  // The numbers above are inconsistent with the equivalent in electron ExcitationEmfietzoglou
  // listed below; this difference must be clarified
  //  const G4double ej[] ={ 8.22*eV, 10.00*eV, 11.24*eV, 12.61*eV, 13.77*eV};

  G4double e = ej[excitationLevel];
  
  return e;
}




//G4bool G4CrossSectionExcitationMillerGreen::IsApplicable(const G4ParticleDefinition& particle)
//{
//  return ( &particle == G4Electron::Electron() );
//} 

