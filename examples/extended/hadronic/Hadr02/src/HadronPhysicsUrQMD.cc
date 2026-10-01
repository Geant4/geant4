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
/// \file HadronPhysicsUrQMD.cc
/// \brief Implementation of the HadronPhysicsUrQMD class

//---------------------------------------------------------------------------
//
// ClassName:
//
// Author: 2012 A. Dotti
//   created from HadronPhysicsUrQMD
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifdef G4_USE_URQMD
#  include "HadronPhysicsUrQMD.hh"

#  include "G4BaryonConstructor.hh"
#  include "G4ChipsKaonMinusInelasticXS.hh"
#  include "G4ChipsKaonPlusInelasticXS.hh"
#  include "G4ChipsKaonZeroInelasticXS.hh"
#  include "G4CrossSectionDataSetRegistry.hh"
#  include "G4VCrossSectionDataSet.hh"
#  include "G4PhysListUtil.hh"
#  include "G4KaonMinus.hh"
#  include "G4KaonPlus.hh"
#  include "G4KaonZeroShort.hh"
#  include "G4KaonZeroLong.hh"
#  include "G4Proton.hh"
#  include "G4Neutron.hh"
#  include "G4PionPlus.hh"
#  include "G4PionMinus.hh"
#  include "G4BGGNucleonInelasticXS.hh"
#  include "G4BGGPionInelasticXS.hh"
#  include "G4HadParticles.hh"
#  include "G4HadronInelasticProcess.hh"
#  include "G4PhysicsListHelper.hh"
#  include "G4MesonConstructor.hh"
#  include "G4ParticleDefinition.hh"
#  include "G4ParticleTable.hh"
#  include "G4PhysicalConstants.hh"
#  include "G4ProcessManager.hh"
#  include "G4ShortLivedConstructor.hh"
#  include "G4SystemOfUnits.hh"
#  include "G4ios.hh"
#  include "globals.hh"

#  include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HadronPhysicsUrQMD::HadronPhysicsUrQMD(G4int) : G4VPhysicsConstructor("hInelastic UrQMD")
{
  fNeutrons = 0;
  fUrQMDNeutron = 0;
  fPiK = 0;
  fUrQMDPiK = 0;
  fPro = 0;
  fUrQMDPro = 0;
  fHyperon = 0;
  fFTFPHyperon = 0;
  fUrQMDAntiBaryon = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HadronPhysicsUrQMD::CreateModels()
{
  fNeutrons = new G4NeutronBuilder;
  fUrQMDNeutron = new UrQMDNeutronBuilder();  // put fission and capture here
  fNeutrons->RegisterMe(fUrQMDNeutron);

  fPro = new G4ProtonBuilder;
  fUrQMDPro = new UrQMDProtonBuilder();
  fPro->RegisterMe(fUrQMDPro);

  fPiK = new G4PiKBuilder;
  fUrQMDPiK = new UrQMDPiKBuilder();
  fPiK->RegisterMe(fUrQMDPiK);

  // For Hyperons use FTF model
  fHyperon = new G4HyperonBuilder;
  fFTFPHyperon = new G4HyperonFTFPBuilder;
  fHyperon->RegisterMe(fFTFPHyperon);

  // Anti-baryons via UrQMD. G4AntiBarionBuilder no longer accepts a custom
  // sub-builder (RegisterMe was removed in 11.x), so the anti-baryon inelastic
  // processes are built directly in ConstructProcess() using this sub-builder.
  fUrQMDAntiBaryon = new UrQMDAntiBarionBuilder();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HadronPhysicsUrQMD::~HadronPhysicsUrQMD()
{
  delete fNeutrons;
  delete fUrQMDNeutron;

  delete fPiK;
  delete fUrQMDPiK;

  delete fPro;
  delete fUrQMDPro;

  delete fHyperon;
  delete fFTFPHyperon;
  delete fUrQMDAntiBaryon;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HadronPhysicsUrQMD::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HadronPhysicsUrQMD::ConstructProcess()
{
  CreateModels();
  fNeutrons->Build();
  fPro->Build();
  fPiK->Build();

  // Modern base builders delegate cross-section registration to the sub-builders,
  // but the UrQMD sub-builders only register the model. Add the standard inelastic
  // cross sections here so the processes have data (fixes had001 "no cross section").
  G4PhysListUtil::FindInelasticProcess(G4Proton::Proton())
    ->AddDataSet(new G4BGGNucleonInelasticXS(G4Proton::Proton()));
  G4PhysListUtil::FindInelasticProcess(G4Neutron::Neutron())
    ->AddDataSet(new G4BGGNucleonInelasticXS(G4Neutron::Neutron()));
  G4PhysListUtil::FindInelasticProcess(G4PionPlus::PionPlus())
    ->AddDataSet(new G4BGGPionInelasticXS(G4PionPlus::PionPlus()));
  G4PhysListUtil::FindInelasticProcess(G4PionMinus::PionMinus())
    ->AddDataSet(new G4BGGPionInelasticXS(G4PionMinus::PionMinus()));

  // use CHIPS cross sections also for Kaons (local datasets)
  G4VCrossSectionDataSet* ChipsKaonMinus =
    G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(
      G4ChipsKaonMinusInelasticXS::Default_Name());
  G4VCrossSectionDataSet* ChipsKaonPlus =
    G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(
      G4ChipsKaonPlusInelasticXS::Default_Name());
  G4VCrossSectionDataSet* ChipsKaonZero =
    G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(
      G4ChipsKaonZeroInelasticXS::Default_Name());
  //

  G4PhysListUtil::FindInelasticProcess(G4KaonMinus::KaonMinus())->AddDataSet(ChipsKaonMinus);
  G4PhysListUtil::FindInelasticProcess(G4KaonPlus::KaonPlus())->AddDataSet(ChipsKaonPlus);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroShort::KaonZeroShort())->AddDataSet(ChipsKaonZero);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroLong::KaonZeroLong())->AddDataSet(ChipsKaonZero);

  fHyperon->Build();

  // Anti-baryons via UrQMD (faithful to the original list; the default
  // G4AntiBarionBuilder would use FTFP). Build each anti-baryon inelastic
  // process directly and register the UrQMD model + anti-nucleon cross section.
  // Set covers pbar, nbar and the light anti-ions (anti d/t/He3/alpha).
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  for (auto pdg : G4HadParticles::GetLightAntiIons()) {
    G4ParticleDefinition* part = table->FindParticle(pdg);
    if (part == nullptr) continue;
    auto* hadi = new G4HadronInelasticProcess(part->GetParticleName() + "Inelastic", part);
    fUrQMDAntiBaryon->Build(hadi);  // registers UrQMD model + anti-nucleon XS
    ph->RegisterProcess(hadi, part);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4HadronicProcess* HadronPhysicsUrQMD::FindInelasticProcess(const G4ParticleDefinition* p)
{
  G4HadronicProcess* had = 0;
  if (p) {
    G4ProcessVector* pvec = p->GetProcessManager()->GetProcessList();
    size_t n = pvec->size();
    if (0 < n) {
      for (size_t i = 0; i < n; ++i) {
        if (fHadronInelastic == ((*pvec)[i])->GetProcessSubType()) {
          had = static_cast<G4HadronicProcess*>((*pvec)[i]);
          break;
        }
      }
    }
  }
  return had;
}

#endif  // G4_USE_URQMD
