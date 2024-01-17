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
// V. Ivanchenko 20 October 2023 
//

#include "G4ParticleHPInelasticXS.hh"
#include "G4ParticleHPManager.hh"
#include "G4SystemOfUnits.hh"

G4ParticleHPInelasticXS::G4ParticleHPInelasticXS(const G4ParticleDefinition* p)
  : G4CrossSectionHP(p, p->GetParticleName() + "InelasticHP",
		     G4ParticleHPManager::GetInstance()->GetParticleHPPath(p) + "/Inelastic/CrossSection/",
                     200*CLHEP::MeV, 0, 100), part(p)
{
  // is the current default, while data should be applicable up to 200 MeV
  SetMaxKinEnergy(30*CLHEP::MeV);
}

void G4ParticleHPInelasticXS::CrossSectionDescription(std::ostream& outF) const
{
  outF << "High Precision cross data based on Evaluated Nuclear Data Files"
       << " (JEFF) for inelastic reaction of " << part->GetParticleName()
       << " below 200 MeV";
}
