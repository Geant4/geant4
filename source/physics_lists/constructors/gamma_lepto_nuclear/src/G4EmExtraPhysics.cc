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
//---------------------------------------------------------------------------
//
// ClassName:   G4EmExtraPhysics
//
// Author: 2002 J.P. Wellisch
//
// Modified:
//
// 10.11.2005 V.Ivanchenko edit to provide a standard
// 19.06.2006 V.Ivanchenko add mu-nuclear process
// 16.10.2012 A.Ribon: renamed G4EmExtraBertiniPhysics as G4EmExtraPhysics
// 10.04.2014 A.Dotti: Add MT functionality for messenger
// 24.04.2014 A.Ribon: switched on muon-nuclear by default
// 29.01.2018 V.Grichine, adding neutrinos 
//
//
///////////////////////////////////////////////////////////////

#include "G4EmExtraPhysics.hh"

#include "G4SystemOfUnits.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoTau.hh"
#include "G4NeutrinoTau.hh"

#include "G4SynchrotronRadiation.hh"
#include "G4BertiniElectroNuclearBuilder.hh"
#include "G4LENDBertiniGammaElectroNuclearBuilder.hh"
#include "G4MuonNuclearProcess.hh"
#include "G4MuonVDNuclearModel.hh"

#include "G4GammaConversionToMuons.hh"
#include "G4AnnihiToMuPair.hh"
#include "G4eeToHadrons.hh"

#include "G4NeutrinoElectronProcess.hh"
#include "G4NeutrinoElectronTotXsc.hh"
#include "G4NeutrinoElectronCcModel.hh"
#include "G4NeutrinoElectronNcModel.hh"
#include "G4GammaGeneralProcess.hh"
#include "G4LossTableManager.hh"

#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"
#include "G4AutoDelete.hh"
 
// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmExtraPhysics);

G4bool G4EmExtraPhysics::gnActivated  = true;
G4bool G4EmExtraPhysics::eActivated   = true;
G4bool G4EmExtraPhysics::gLENDActivated = false;
G4bool G4EmExtraPhysics::munActivated = true;
G4bool G4EmExtraPhysics::synActivated = false;
G4bool G4EmExtraPhysics::synActivatedForAll = false;
G4bool G4EmExtraPhysics::gmumuActivated = false;
G4bool G4EmExtraPhysics::pmumuActivated = false;
G4bool G4EmExtraPhysics::phadActivated  = false;
G4bool G4EmExtraPhysics::fNuActivated  = false;
G4bool G4EmExtraPhysics::fNuETotXscActivated  = false;

G4double G4EmExtraPhysics::gmumuFactor  = 1.0;
G4double G4EmExtraPhysics::pmumuFactor  = 1.0;
G4double G4EmExtraPhysics::phadFactor   = 1.0;
G4double G4EmExtraPhysics::fNuEleCcBias = 1.0;
G4double G4EmExtraPhysics::fNuEleNcBias = 1.0;
G4double G4EmExtraPhysics::fNuNucleusBias = 1.0;

G4String G4EmExtraPhysics::fNuDetectorName = "0";

G4ThreadLocal G4BertiniElectroNuclearBuilder* G4EmExtraPhysics::theGNPhysics = nullptr;
G4ThreadLocal G4SynchrotronRadiation* G4EmExtraPhysics::theSynchRad = nullptr;
G4ThreadLocal G4GammaConversionToMuons* G4EmExtraPhysics::theGammaToMuMu = nullptr;
G4ThreadLocal G4AnnihiToMuPair* G4EmExtraPhysics::thePosiToMuMu = nullptr;
G4ThreadLocal G4eeToHadrons* G4EmExtraPhysics::thePosiToHadrons = nullptr;

G4ThreadLocal G4NeutrinoElectronProcess* G4EmExtraPhysics::theNuEleProcess = nullptr;
G4ThreadLocal G4NeutrinoElectronTotXsc*  G4EmExtraPhysics::theNuEleTotXsc  = nullptr;


//////////////////////////////////////

G4EmExtraPhysics::G4EmExtraPhysics(G4int ver): 
  G4VPhysicsConstructor("G4GammaLeptoNuclearPhys"),
  verbose(ver)
{
  theMessenger = new G4EmMessenger(this);
  SetPhysicsType(bEmExtra);
  if(verbose > 1) G4cout << "### G4EmExtraPhysics" << G4endl;
}

G4EmExtraPhysics::G4EmExtraPhysics(const G4String&)
  : G4EmExtraPhysics(1)
{}

G4EmExtraPhysics::~G4EmExtraPhysics()
{
  delete theMessenger;
  theMessenger = nullptr;
}

void G4EmExtraPhysics::Synch(G4bool val)
{
  synActivated = val;
}

void G4EmExtraPhysics::SynchAll(G4bool val)
{
  synActivatedForAll = val;
  if(synActivatedForAll) { synActivated = true; }
}

void G4EmExtraPhysics::GammaNuclear(G4bool val)
{
  gnActivated = val;
}

void G4EmExtraPhysics::LENDGammaNuclear(G4bool val)
{
  gLENDActivated = val;
}

void G4EmExtraPhysics::ElectroNuclear(G4bool val)
{
  eActivated = val;
}

void G4EmExtraPhysics::MuonNuclear(G4bool val)
{
  munActivated = val;
}

void G4EmExtraPhysics::GammaToMuMu(G4bool val)
{
  gmumuActivated = val;
}

void G4EmExtraPhysics::PositronToMuMu(G4bool val)
{
  pmumuActivated = val;
}

void G4EmExtraPhysics::PositronToHadrons(G4bool val)
{
  phadActivated = val;
}

void G4EmExtraPhysics::GammaToMuMuFactor(G4double val)
{
  if(val > 0.0) gmumuFactor = val;
}

void G4EmExtraPhysics::PositronToMuMuFactor(G4double val)
{
  if(val > 0.0) pmumuFactor = val;
}

void G4EmExtraPhysics::PositronToHadronsFactor(G4double val)
{
  if(val > 0.0) phadFactor = val;
}

////////////////////////////////////////////////////

void G4EmExtraPhysics::NeutrinoActivated(G4bool val)
{
  fNuActivated = val;
}

void G4EmExtraPhysics::NuETotXscActivated(G4bool val)
{
  fNuETotXscActivated = val;
}

void G4EmExtraPhysics::SetNuEleCcBias(G4double bf)
{
  if(bf > 0.0) fNuEleCcBias = bf;
}

void G4EmExtraPhysics::SetNuEleNcBias(G4double bf)
{
  if(bf > 0.0) fNuEleNcBias = bf;
}

void G4EmExtraPhysics::SetNuNucleusBias(G4double bf)
{
  if(bf > 0.0) fNuNucleusBias = bf;
}

void G4EmExtraPhysics::SetNuDetectorName(const G4String& dn)
{
  fNuDetectorName = dn;
}

/////////////////////////////////////////////////

void G4EmExtraPhysics::ConstructParticle()
{
  G4Gamma::Gamma();
  G4Electron::Electron();
  G4Positron::Positron();
  G4MuonPlus::MuonPlus();
  G4MuonMinus::MuonMinus();

  G4AntiNeutrinoE::AntiNeutrinoE();
  G4NeutrinoE::NeutrinoE();
  G4AntiNeutrinoMu::AntiNeutrinoMu();
  G4NeutrinoMu::NeutrinoMu();
  G4AntiNeutrinoTau::AntiNeutrinoTau();
  G4NeutrinoTau::NeutrinoTau();
}

void G4EmExtraPhysics::ConstructProcess()
{
  G4ParticleDefinition* gamma     = G4Gamma::Gamma();
  G4ParticleDefinition* electron  = G4Electron::Electron();
  G4ParticleDefinition* positron  = G4Positron::Positron();
  G4ParticleDefinition* muonplus  = G4MuonPlus::MuonPlus();
  G4ParticleDefinition* muonminus = G4MuonMinus::MuonMinus();

  G4ParticleDefinition* anuelectron = G4AntiNeutrinoE::AntiNeutrinoE();
  G4ParticleDefinition* nuelectron  = G4NeutrinoE::NeutrinoE(); 
  G4ParticleDefinition* anumuon     = G4AntiNeutrinoMu::AntiNeutrinoMu();
  G4ParticleDefinition* numuon      = G4NeutrinoMu::NeutrinoMu();
  G4ParticleDefinition* anutau      = G4AntiNeutrinoTau::AntiNeutrinoTau();
  G4ParticleDefinition* nutau       = G4NeutrinoTau::NeutrinoTau();

  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  if(gnActivated) {
    if ( gLENDActivated != true ) {
      theGNPhysics = new G4BertiniElectroNuclearBuilder(eActivated);
    } else {
      theGNPhysics = new G4LENDBertiniGammaElectroNuclearBuilder(eActivated);
    }
    theGNPhysics->Build();
  }
  if(munActivated) {
    G4MuonNuclearProcess* muNucProcess = new G4MuonNuclearProcess();
    G4MuonVDNuclearModel* muNucModel = new G4MuonVDNuclearModel();
    muNucProcess->RegisterMe(muNucModel);
    ph->RegisterProcess( muNucProcess, muonplus);
    ph->RegisterProcess( muNucProcess, muonminus);
  }
  if(gmumuActivated) {
    theGammaToMuMu = new G4GammaConversionToMuons();
    theGammaToMuMu->SetCrossSecFactor(gmumuFactor);
    G4GammaGeneralProcess* sp = 
      (G4GammaGeneralProcess*)G4LossTableManager::Instance()->GetGammaGeneralProcess();
    if(sp) {
      sp->AddMMProcess(theGammaToMuMu);
    } else {
      ph->RegisterProcess(theGammaToMuMu, gamma);
    }
  }  
  if(pmumuActivated) {
    thePosiToMuMu = new G4AnnihiToMuPair();
    thePosiToMuMu->SetCrossSecFactor(pmumuFactor);
    ph->RegisterProcess(thePosiToMuMu, positron);
  }  
  if(phadActivated) {
    thePosiToHadrons = new G4eeToHadrons();
    thePosiToHadrons->SetCrossSecFactor(phadFactor);
    ph->RegisterProcess(thePosiToHadrons, positron);
  }  
  if(synActivated) {
    theSynchRad = new G4SynchrotronRadiation();
    ph->RegisterProcess( theSynchRad, electron);
    ph->RegisterProcess( theSynchRad, positron);
    if(synActivatedForAll) {
      auto myParticleIterator=GetParticleIterator();
      myParticleIterator->reset();
      G4ParticleDefinition* particle = nullptr;

      while( (*myParticleIterator)() ) {
	particle = myParticleIterator->value();
	if( particle->GetPDGStable() && particle->GetPDGCharge() != 0.0) { 
	  if(verbose > 1) {
	    G4cout << "### G4SynchrotronRadiation for " 
		   << particle->GetParticleName() << G4endl;
	  }
	  ph->RegisterProcess( theSynchRad, particle);
	}
      }
    }
  }
  if( fNuActivated )
  {
    theNuEleProcess = new G4NeutrinoElectronProcess(fNuDetectorName);
    theNuEleTotXsc = new G4NeutrinoElectronTotXsc();

    if(fNuETotXscActivated)
    {
      G4double bftot = std::max(fNuEleCcBias,fNuEleNcBias);
      theNuEleProcess->SetBiasingFactor(bftot);
    }
    else
    {
      theNuEleProcess->SetBiasingFactors(fNuEleCcBias,fNuEleNcBias);
      theNuEleTotXsc->SetBiasingFactors(fNuEleCcBias,fNuEleNcBias);
    }
    theNuEleProcess->AddDataSet(theNuEleTotXsc);

    G4NeutrinoElectronCcModel* ccModel = new G4NeutrinoElectronCcModel();
    G4NeutrinoElectronNcModel* ncModel = new G4NeutrinoElectronNcModel();
    theNuEleProcess->RegisterMe(ccModel);
    theNuEleProcess->RegisterMe(ncModel);

    ph->RegisterProcess(theNuEleProcess, anuelectron);
    ph->RegisterProcess(theNuEleProcess, nuelectron);
    ph->RegisterProcess(theNuEleProcess, anumuon);
    ph->RegisterProcess(theNuEleProcess, numuon);
    ph->RegisterProcess(theNuEleProcess, anutau);
    ph->RegisterProcess(theNuEleProcess, nutau);
  }
}

