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
//
//
// G4 Physics class: HadronElasticDataSet for cross sections
// F.W. Jones, TRIUMF, 28-JAN-98
// 
// Modified: V.Ivanchenko
//

#include "G4HadronElasticDataSet.hh"
#include "G4HadronCrossSections.hh"
#include "G4DynamicParticle.hh"
#include "G4NistManager.hh"
#include <iostream>

G4HadronElasticDataSet::G4HadronElasticDataSet(const G4String& nam)
  : G4VCrossSectionDataSet(nam), theZ(0),fElasticXS(0.0),
    fKinEnergy(0.0),fParticle(nullptr)
{
  fGheishaXS = G4HadronCrossSections::Instance(); 
  fNIST = G4NistManager::Instance();
}

G4HadronElasticDataSet::~G4HadronElasticDataSet() {}

void 
G4HadronElasticDataSet::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4HadronElasticDataSet contains elastic cross sections for\n"
          << "all long-lived hadrons at all incident energies.  It was\n"
          << "developed as part of the Gheisha hadronic package\n"
          << "by H. Fesefeldt, and consists of a set of parameterizations\n"
          << "of elastic scattering data.\n";
}

G4bool
G4HadronElasticDataSet::IsElementApplicable(const G4DynamicParticle*, 
					    G4int, const G4Material*)
{
  return true;
}

G4double G4HadronElasticDataSet::GetElementCrossSection(
         const G4DynamicParticle* aParticle, G4int Z, const G4Material*)
{
  G4double ekin = aParticle->GetKineticEnergy();
  const G4ParticleDefinition* pd = aParticle->GetDefinition();
  if(Z != theZ || ekin != fKinEnergy || pd != fParticle) {
    theZ = Z;
    fKinEnergy = ekin;
    fParticle = pd;
    G4int A = G4lrint(fNIST->GetAtomicMassAmu(Z));
    fElasticXS = fGheishaXS->GetElasticCrossSection(aParticle, Z, A);
  }
  return fElasticXS;
}
