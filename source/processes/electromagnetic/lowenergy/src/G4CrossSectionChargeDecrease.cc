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
// $Id: G4CrossSectionChargeDecrease.cc,v 1.7 2009-06-11 15:47:08 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4CrossSectionChargeDecrease.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CrossSectionChargeDecrease::G4CrossSectionChargeDecrease()
{
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
    highEnergyLimit[proton] = 10. * MeV;
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

   G4cout << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "   The class G4CrossSectionChargeDecrease is NOT SUPPORTED ANYMORE. " << G4endl;
   G4cout << "   It will be REMOVED with the next major release of Geant4. " << G4endl;
   G4cout << "   Please consult: https://twiki.cern.ch/twiki/bin/view/Geant4/LoweProcesses" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CrossSectionChargeDecrease::~G4CrossSectionChargeDecrease()
{}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

  std::map< G4String,G4double,std::less<G4String> >::iterator pos1;
  pos1 = lowEnergyLimit.find(particleName);

  if (pos1 != lowEnergyLimit.end())
  {
    lowLim = pos1->second;
  }

  std::map< G4String,G4double,std::less<G4String> >::iterator pos2;
  pos2 = highEnergyLimit.find(particleName);

  if (pos2 != highEnergyLimit.end())
  {
    highLim = pos2->second;
  }

  G4double crossSection(0.);
  if (k >= lowLim && k <= highLim)
  {
      crossSection = partialCrossSection.Sum(k,particleDefinition);
  }
 
  return crossSection;
}

