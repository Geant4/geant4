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
// $Id: G4CrossSectionChargeIncrease.cc,v 1.7 2009-06-11 15:47:08 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4CrossSectionChargeIncrease.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CrossSectionChargeIncrease::G4CrossSectionChargeIncrease()
{
  lowEnergyLimitDefault = 1 * keV;
  highEnergyLimitDefault = 10 * MeV;

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();
  G4ParticleDefinition* hydrogenDef = instance->GetIon("hydrogen");
  G4ParticleDefinition* alphaPlusDef = instance->GetIon("alpha+");
  G4ParticleDefinition* heliumDef = instance->GetIon("helium");

  G4String hydrogen;
  G4String alphaPlus;
  G4String helium;

  if (hydrogenDef != 0)
  {
    hydrogen = hydrogenDef->GetParticleName();
    lowEnergyLimit[hydrogen] = 1. * keV;
    highEnergyLimit[hydrogen] = 10. * MeV;
  }
  else
  {
    G4Exception("G4CrossSectionChargeIncrease Constructor: hydrogen is not defined");
  }

  if (alphaPlusDef != 0)
  {
    alphaPlus = alphaPlusDef->GetParticleName();
    lowEnergyLimit[alphaPlus] = 1. * keV;
    highEnergyLimit[alphaPlus] = 10. * MeV;
  }
  else
  {
    G4Exception("G4CrossSectionChargeIncrease Constructor: alphaPlus is not defined");
  }

  if (heliumDef != 0)
  {
    helium = heliumDef->GetParticleName();
    lowEnergyLimit[helium] = 1. * keV;
    highEnergyLimit[helium] = 10. * MeV;
  }
  else
  {
    G4Exception("G4CrossSectionChargeIncrease Constructor: helium is not defined");
  }

   G4cout << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "   The class G4CrossSectionChargeIncrease is NOT SUPPORTED ANYMORE. " << G4endl;
   G4cout << "   It will be REMOVED with the next major release of Geant4. " << G4endl;
   G4cout << "   Please consult: https://twiki.cern.ch/twiki/bin/view/Geant4/LoweProcesses" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CrossSectionChargeIncrease::~G4CrossSectionChargeIncrease()
{}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4CrossSectionChargeIncrease::CrossSection(const G4Track& track)
{
  G4double lowLim = lowEnergyLimitDefault;
  G4double highLim = highEnergyLimitDefault;

  const G4DynamicParticle* particle = track.GetDynamicParticle();
  G4double k = particle->GetKineticEnergy();

  const G4String& particleName = particle->GetDefinition()->GetParticleName();

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();
  
  const G4ParticleDefinition* particleDefinition = track.GetDefinition();

  if (
      particleDefinition != instance->GetIon("hydrogen")
      &&
      particleDefinition != instance->GetIon("alpha+")
      &&
      particleDefinition != instance->GetIon("helium")
      )
   	    
    G4Exception("G4CrossSectionChargeIncrease: attempting to calculate cross section for wrong particle");

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

  G4double totalCrossSection = 0.;

  if (k >= lowLim && k <= highLim)
  {
      //HYDROGEN
      if (particleDefinition == instance->GetIon("hydrogen"))
      {
	  const  G4double aa = 2.835;
	  const  G4double bb = 0.310;
	  const  G4double cc = 2.100;
	  const  G4double dd = 0.760;
	  const  G4double fac = 1.0e-18;
	  const  G4double rr = 13.606 * eV;  
	  
	  G4double t = k / (proton_mass_c2/electron_mass_c2); 
	  G4double x = t / rr; 
	  G4double temp = 4.0 * pi * Bohr_radius/nm * Bohr_radius/nm * fac;
	  G4double sigmal =  temp * cc * (std::pow(x,dd));
	  G4double sigmah = temp * (aa * std::log(1.0 + x) + bb) / x;
	  totalCrossSection = 1.0/(1.0/sigmal + 1.0/sigmah) *m*m;
      }
      else
      {
	  totalCrossSection = partialCrossSection.Sum(k,particleDefinition);
      }
  }
  
  return totalCrossSection;
}

