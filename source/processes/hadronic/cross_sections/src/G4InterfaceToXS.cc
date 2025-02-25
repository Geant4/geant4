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
// GEANT4 Class file
//
// File name:  G4InterfaceToXS
//
// Authors  V.Ivantchenko  12 July 2024
//
// The interface provide access to inelastic cross section of gamma,
// neutrons, and light ions based on G4PARTICLEXS evaluated data.
// The only one instance of the cross section class exist in each thread.
// This class may be used by any model, cross section or other classes.

#include "G4InterfaceToXS.hh"
#include "G4ParticleDefinition.hh"
#include "G4GammaNuclearXS.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4ParticleInelasticXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4Log.hh"

G4InterfaceToXS::G4InterfaceToXS(const G4ParticleDefinition* p, G4int idx)
  : index(idx), fParticle(p)
{
  auto reg = G4CrossSectionDataSetRegistry::Instance();
  if (index == 0) {
    fNeutronNuclear =
      dynamic_cast<G4NeutronInelasticXS*>(reg->GetCrossSectionDataSet("G4NeutronInelasticXS"));
    if (nullptr == fNeutronNuclear) { fNeutronNuclear = new G4NeutronInelasticXS(); }
    fNeutronNuclear->BuildPhysicsTable(*fParticle);
  } else if (index > 0 && index < 6) {
    const G4String pname[6] = {"neutron", "proton", "deuteron", "triton", "He3", "alpha"};
    const G4String& ss = pname[index] + "ParticleXS";
    fParticleNuclear =
      dynamic_cast<G4ParticleInelasticXS*>(reg->GetCrossSectionDataSet(ss));
    if (nullptr == fParticleNuclear) { fParticleNuclear = new G4ParticleInelasticXS(fParticle); }
    fParticleNuclear->BuildPhysicsTable(*fParticle);
  } else if (index == 6) {
    fGammaNuclear =
      dynamic_cast<G4GammaNuclearXS*>(reg->GetCrossSectionDataSet("GammaNuclearXS"));
    if (nullptr == fGammaNuclear) { fGammaNuclear = new G4GammaNuclearXS(); }
    fGammaNuclear->BuildPhysicsTable(*fParticle);
  }
}

G4double G4InterfaceToXS::GetElementCrossSection(const G4double ekin, const G4int Z) 
{
  G4double res = 0.0;
  if (ekin <= 0.0) { return res; }
  if (nullptr != fNeutronNuclear) {
    res = fNeutronNuclear->ElementCrossSection(ekin, G4Log(ekin), Z);
  } else if (nullptr != fParticleNuclear) {
    res = fParticleNuclear->ElementCrossSection(ekin, G4Log(ekin), Z);
  } else if (nullptr != fGammaNuclear) {
    res = fGammaNuclear->ElementCrossSection(ekin, Z);
  }
  return res;
}

G4double G4InterfaceToXS::GetIsoCrossSection(const G4double ekin, const G4int Z, const G4int A) 
{
  G4double res = 0.0;
  if (ekin <= 0.0) { return res; }
  if (nullptr != fNeutronNuclear) {
    res = fNeutronNuclear->IsoCrossSection(ekin, G4Log(ekin), Z, A);
  } else if (nullptr != fParticleNuclear) {
    res = fParticleNuclear->IsoCrossSection(ekin, G4Log(ekin), Z, A);
  } else if (nullptr != fGammaNuclear) {
    res = fGammaNuclear->IsoCrossSection(ekin, Z, A);
  }
  return res;
}
