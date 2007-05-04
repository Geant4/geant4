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
// $Id: G4eCrossSectionExcitationEmfietzoglou.cc,v 1.1 2007-05-04 10:16:06 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// $Id: G4eCrossSectionExcitationEmfietzoglou.cc,v 1.1 2007-05-04 10:16:06 pia Exp $
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


#include "G4eCrossSectionExcitationEmfietzoglou.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"

#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Electron.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"

#include "Randomize.hh"

G4eCrossSectionExcitationEmfietzoglou::G4eCrossSectionExcitationEmfietzoglou()
{

  name = "eCrossSectionExcitationEmfietzoglou";
  lowEnergyLimit = 7. * eV;
  highEnergyLimit = 10 * keV;

//  if (verboseLevel > 0)
//  {
//    G4cout << name << " is created " << G4endl
//     << "Energy range: "
//     << lowEnergyLimit / keV << " keV - "
//     << highEnergyLimit / GeV << " GeV"
//     << G4endl;
//  }
}


G4eCrossSectionExcitationEmfietzoglou::~G4eCrossSectionExcitationEmfietzoglou()
{ }
 

G4double G4eCrossSectionExcitationEmfietzoglou::CrossSection(const G4Track& track)
{
  const G4DynamicParticle* particle = track.GetDynamicParticle();
  G4double k = particle->GetKineticEnergy();

  // G4Material* material = track.GetMaterial();

  // Assume that the material is water; proper algorithm to calculate correctly for any material to be inserted here

  // Take into account 5 excitation levels (D. Emfietzoglou et al., NIM B 193, pp. 71-78, 2002.
  G4int i = 5;

  G4double crossSection(0.);

  while (i>0)
    {
      i--;
      crossSection += PartialCrossSection(k,i);
    }
  return crossSection;
 }


G4double G4eCrossSectionExcitationEmfietzoglou::PartialCrossSection(G4double k, G4int excitationLevel)
{
  //                 Aj                        T
  // sigma(T) = ------------- (Bj /  T) ln(Cj ---) [1 - Bj / T]^Pj
  //             2 pi alpha0                   R
  //
  // T      is the incoming electron kinetic energy
  // alpha0 is the Bohr Radius (Bohr_radius)
  // Aj, Bj, Cj & Pj are parameters that can be found in Emfietzoglou's papers
  //
  //
  // From Phys. Med. Biol. 48 (2003) 2355-2371, D.Emfietzoglou,
  // Monte Carlo Simulation of the energy loss of low energy electrons in liquid Water
  
  k = k / eV;
  const G4double sigma0 = (10. / 3.343e22) * cm2;
  const G4double aj[] = {0.0205, 0.0209, 0.0130, 0.0026, 0.0025};
  const G4double cj[] = {4.9801, 3.3850, 2.8095, 1.9242, 3.4624};
  const G4double pj[] = {0.4757, 0.3483, 0.4443, 0.3429, 0.4379};
  const G4double r = 13.6 * eV;
  
  G4double cross = 0.;
  if (k < EnergyConstant(excitationLevel)) return cross; 

  G4double excitationSigma = ( aj[excitationLevel] / (2. * pi * Bohr_radius)) *
    (EnergyConstant(excitationLevel) / k) * 
    std::log(cj[excitationLevel] * (k / r)) * 
    std::pow((1. - (EnergyConstant(excitationLevel) / k)), pj[excitationLevel]);
  G4cout << "Bohr radius = " << Bohr_radius << G4endl;

  cross = excitationSigma * sigma0;

  return cross;
}


G4int G4eCrossSectionExcitationEmfietzoglou::RandomizePartialCrossSection(G4double k)
{
  // Assume 5 excitation levels 

  G4int i = 5;
  G4double value = 0.;
  G4double values[5];
  
  while (i>0)
    {
      i--;
      values[i] = PartialCrossSection(k,i);
      value += values[i];
    }
  
  value *= G4UniformRand();
  
  i = 5;
  while (i>0)
    {
      i--;
      
      if (values[i] > value) return i;
      
      value -= values[i];
    }
  // One should never end up here; next statement added to avoid compilation warning
  // Probably one should throw an exception if one ends up here
  return 0;
}

G4double G4eCrossSectionExcitationEmfietzoglou::EnergyConstant(G4int excitationLevel)
{
  const G4double ej[] ={ 8.22*eV, 10.00*eV, 11.24*eV, 12.61*eV, 13.77*eV};
  G4double e = ej[excitationLevel];
  
  return e;
}




//G4bool G4eCrossSectionExcitationEmfietzoglou::IsApplicable(const G4ParticleDefinition& particle)
//{
//  return ( &particle == G4Electron::Electron() );
//} 

