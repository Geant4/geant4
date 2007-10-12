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
// $Id: G4eCrossSectionScreenedRutherford.cc,v 1.3 2007-10-12 12:27:19 pia Exp $
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


#include "G4eCrossSectionScreenedRutherford.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"


G4eCrossSectionScreenedRutherford::G4eCrossSectionScreenedRutherford()
{

  name = "eCrossSectionScreenedRutherford";
  lowEnergyLimit = 7. * eV;
  highEnergyLimit = 10 * MeV;

//  if (verboseLevel > 0)
//  {
//    G4cout << name << " is created " << G4endl
//     << "Energy range: "
//     << lowEnergyLimit / keV << " keV - "
//     << highEnergyLimit / GeV << " GeV"
//     << G4endl;
//  }
}


G4eCrossSectionScreenedRutherford::~G4eCrossSectionScreenedRutherford()
{ }
 

G4double G4eCrossSectionScreenedRutherford::CrossSection(const G4Track& track)
{
  const G4DynamicParticle* particle = track.GetDynamicParticle();
  G4double k = particle->GetKineticEnergy();

  // Cross section = 0 outside the energy validity limits set in the constructor
  // ---- MGP ---- Better handling of these limits to be set in a following design iteration 

  G4double screenedCrossSection = 0.;

  if (k > lowEnergyLimit && k < highEnergyLimit)
    {      
      // G4Material* material = track.GetMaterial();

      // Assume that the material is water; proper algorithm to calculate z correctly for any material to be inserted here
      // For H20 Z = 10 (total number of electrons)
      G4double z = 10.;
      
      G4double n = ScreeningFactor(k,z);
      G4double crossSection = RutherfordCrossSection(k, z);
      screenedCrossSection = pi *  crossSection / (n * (n + 1.));
    }    

  return screenedCrossSection;
}
  
G4double G4eCrossSectionScreenedRutherford::RutherfordCrossSection(G4double k, G4double z)
{
  //   
  //                               e^4         /      K + m_e c^2      \^2
  // sigma_Ruth(K) = Z (Z+1) -------------------- | --------------------- |
  //                          (4 pi epsilon_0)^2  \  K * (K + 2 m_e c^2)  /
  //
  // Where K is the electron non-relativistic kinetic energy
  // 
  // NIM 155, pp. 145-156, 1978
  
  G4double length =(e_squared * (k + electron_mass_c2)) / (4 * pi *epsilon0 * k * ( k + 2 * electron_mass_c2));
  G4double cross = z * ( z + 1) * length * length;
  
  return cross;
}

//G4bool G4eCrossSectionScreenedRutherford::IsApplicable(const G4ParticleDefinition& particle)
//{
//  return ( &particle == G4Electron::Electron() );
//} 


G4double G4eCrossSectionScreenedRutherford::ScreeningFactor(G4double k, G4double z)
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
  // NIM 155, pp. 145-156, 1978
  // Formulae (2) and (5)

  const G4double alpha_1(1.64);
  const G4double beta_1(-0.0825);
  const G4double constK(1.7E-5);

  G4double numerator = (alpha_1 + beta_1 * std::log(k/eV)) * constK * std::pow(z, 2./3.);

  k /= electron_mass_c2;

  G4double denominator = k * (2 + k);

  G4double value = 0.;
  if (denominator > 0.) value = numerator / denominator;

  return value;

}

