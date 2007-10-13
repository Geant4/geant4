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
// $Id: G4CrossSectionExcitationEmfietzoglou.cc,v 1.1 2007-10-13 01:55:59 pia Exp $
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
// Reference: TNS Geant4-DNA paper
// S. Chauvie et al., Geant4 physics processes for microdosimetry simulation:
// design foundation and implementation of the first set of models,
// IEEE Trans. Nucl. Sci., vol. 54, no. 6, Dec. 2007.
// Further documentation available from http://www.ge.infn.it/geant4/dna

// -------------------------------------------------------------------


#include "G4CrossSectionExcitationEmfietzoglou.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"


G4CrossSectionExcitationEmfietzoglou::G4CrossSectionExcitationEmfietzoglou()
{

  name = "CrossSectionExcitationEmfietzoglou";
  lowEnergyLimit = 7.4 * eV;
  highEnergyLimit = 10. * MeV;
  nLevels = 5;

  energyConstant.push_back(8.22*eV);
  energyConstant.push_back(10.00*eV);
  energyConstant.push_back(11.24*eV);
  energyConstant.push_back(12.61*eV);
  energyConstant.push_back(13.77*eV);

//  if (verboseLevel > 0)
//  {
//    G4cout << name << " is created " << G4endl
//     << "Energy range: "
//     << lowEnergyLimit / keV << " keV - "
//     << highEnergyLimit / GeV << " GeV"
//     << G4endl;
//  }
}


G4CrossSectionExcitationEmfietzoglou::~G4CrossSectionExcitationEmfietzoglou()
{ }
 

G4double G4CrossSectionExcitationEmfietzoglou::CrossSection(const G4Track& track)
{
  const G4DynamicParticle* particle = track.GetDynamicParticle();
  G4double k = particle->GetKineticEnergy();

  // Cross section = 0 outside the energy validity limits set in the constructor
  // ---- MGP ---- Better handling of these limits to be set in a following design iteration 

  G4double totalCrossSection = 0.;

  if (k > lowEnergyLimit && k < highEnergyLimit)
    {      
      for (G4int i=0; i<nLevels; i++)
	{
	  totalCrossSection += PartialCrossSection(k,i);
	}
    }

  return totalCrossSection;
}

G4double G4CrossSectionExcitationEmfietzoglou::PartialCrossSection(G4double t, G4int level)
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
  
  const G4double sigma0 = (10. / 3.343e22) * cm2;
  
  const G4double aj[]={0.0205, 0.0209, 0.0130, 0.0026, 0.0025};
  const G4double cj[]={4.9801, 3.3850, 2.8095, 1.9242, 3.4624};
  const G4double pj[]={0.4757, 0.3483, 0.4443, 0.3429, 0.4379};
  const G4double r = 13.6 * eV;
  
  G4double sigma = 0.;
  
  if (t >= energyConstant[level])
    {
      G4double excSigma = ( aj[level] / (2.*pi*Bohr_radius)) 
	* (energyConstant[level] / t) 
	* log(cj[level]*(t/r)) 
	* pow((1.- (energyConstant[level]/t)), pj[level]);
      sigma = excSigma * sigma0;
    }
  return sigma;
}



