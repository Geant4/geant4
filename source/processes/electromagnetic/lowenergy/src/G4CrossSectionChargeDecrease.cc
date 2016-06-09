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
// $Id: G4CrossSectionChargeDecrease.cc,v 1.2 2007/11/09 20:11:04 pia Exp $
// GEANT4 tag $Name: geant4-09-01 $
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


#include "G4CrossSectionChargeDecrease.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include "G4Proton.hh"
#include "G4CrossSectionExcitationEmfietzoglouPartial.hh"
#include "G4DNAGenericIonsManager.hh"

G4CrossSectionChargeDecrease::G4CrossSectionChargeDecrease()
{
  // Default energy limits (defined for protection against anomalous behaviour only)
  name = "ChargeDecrease";
  lowEnergyLimitDefault = 1 * keV;
  highEnergyLimitDefault = 10 * MeV;

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();
  G4ParticleDefinition* protonDef = G4Proton::ProtonDefinition();
  G4ParticleDefinition* alphaPlusPlusDef = instance->GetIon("alpha++");
  G4ParticleDefinition* alphaPlusDef = instance->GetIon("alpha+");

  G4String proton;
  G4String alphaPlusPlus;
  G4String alphaPlus;

  if (protonDef != 0)
    {
      proton = protonDef->GetParticleName();
      lowEnergyLimit[proton] = 1. * keV;
      highEnergyLimit[proton] = 10. * keV;
    }
  else
    {
      G4Exception("G4CrossSectionChargeDecrease Constructor: proton is not defined");
    }

  if (alphaPlusPlusDef != 0)
    {
      alphaPlusPlus = alphaPlusPlusDef->GetParticleName();
      lowEnergyLimit[alphaPlusPlus] = 1. * keV;
      highEnergyLimit[alphaPlusPlus] = 10. * MeV;
    }
  else
    {
      G4Exception("G4CrossSectionChargeDecrease Constructor: alphaPlusPlus is not defined");
    }

  if (alphaPlusDef != 0)
    {
      alphaPlus = alphaPlusDef->GetParticleName();
      lowEnergyLimit[alphaPlus] = 1. * keV;
      highEnergyLimit[alphaPlus] = 10. * MeV;
    }
  else
    {
      G4Exception("G4CrossSectionChargeDecrease Constructor: alphaPlus is not defined");
    }

}


G4CrossSectionChargeDecrease::~G4CrossSectionChargeDecrease()
{}
 

G4double G4CrossSectionChargeDecrease::CrossSection(const G4Track& track)
{
  G4double lowLim = lowEnergyLimitDefault;
  G4double highLim = highEnergyLimitDefault;

  const G4DynamicParticle* particle = track.GetDynamicParticle();
  G4double k = particle->GetKineticEnergy();

  const G4ParticleDefinition* particleDefinition = track.GetDefinition();
  
  const G4String& particleName = particleDefinition->GetParticleName();

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  if (
      particleDefinition != G4Proton::ProtonDefinition()
      &&
      particleDefinition != instance->GetIon("alpha++")
      &&
      particleDefinition != instance->GetIon("alpha+")
      )
   	    
    G4Exception("G4CrossSectionChargeDecrease: attempting to calculate cross section for wrong particle");

  // Retrieve energy limits for the current particle type

  std::map< G4String,G4double,std::less<G4String> >::iterator pos1;
  pos1 = lowEnergyLimit.find(particleName);

  // Lower limit
  if (pos1 != lowEnergyLimit.end())
    {
      lowLim = pos1->second;
    }

  // Upper limit
  std::map< G4String,G4double,std::less<G4String> >::iterator pos2;
  pos2 = highEnergyLimit.find(particleName);

  if (pos2 != highEnergyLimit.end())
    {
      highLim = pos2->second;
    }

  //

  G4double crossSection(0.);
  if (k >= lowLim && k <= highLim)
    {
      crossSection = partialCrossSection.Sum(k,particleDefinition);
    }
 
  return crossSection;
}

