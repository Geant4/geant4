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
// 07.05.2019 V.Grichine, adding muon neutrino nucleus interactions 
// 03.11.2022 V. Grichne update for tau-neutrino nucleus processes
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

#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4GenericIon.hh"

#include "G4SynchrotronRadiation.hh"
#include "G4MuonNuclearProcess.hh"
#include "G4MuonVDNuclearModel.hh"
#include "G4ElectroVDNuclearModel.hh"
#include "G4TheoFSGenerator.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4QGSModel.hh"
#include "G4GammaParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4CascadeInterface.hh"
#include "G4LowEGammaNuclearModel.hh"

#include "G4LENDorBERTModel.hh"
#include "G4LENDCombinedCrossSection.hh"

#include "G4GammaConversionToMuons.hh"
#include "G4AnnihiToMuPair.hh"
#include "G4eeToHadrons.hh"
#include "G4MuonToMuonPairProduction.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4ElectronNuclearProcess.hh"
#include "G4PositronNuclearProcess.hh"

#include "G4NeutrinoElectronProcess.hh"
#include "G4NeutrinoElectronTotXsc.hh"
#include "G4NeutrinoElectronCcModel.hh"
#include "G4NeutrinoElectronNcModel.hh"

#include "G4MuNeutrinoNucleusProcess.hh"
#include "G4TauNeutrinoNucleusProcess.hh"
#include "G4ElNeutrinoNucleusProcess.hh"

#include "G4MuNeutrinoNucleusTotXsc.hh"
#include "G4TauNeutrinoNucleusTotXsc.hh"
#include "G4ElNeutrinoNucleusTotXsc.hh"

#include "G4NuMuNucleusCcModel.hh"
#include "G4NuMuNucleusNcModel.hh"
#include "G4ANuMuNucleusCcModel.hh"
#include "G4ANuMuNucleusNcModel.hh"

#include "G4NuTauNucleusCcModel.hh"
#include "G4NuTauNucleusNcModel.hh"
#include "G4ANuTauNucleusCcModel.hh"
#include "G4ANuTauNucleusNcModel.hh"

#include "G4NuElNucleusCcModel.hh"
#include "G4NuElNucleusNcModel.hh"
#include "G4ANuElNucleusCcModel.hh"
#include "G4ANuElNucleusNcModel.hh"

#include "G4GammaGeneralProcess.hh"
#include "G4LossTableManager.hh"
#include "G4PhotoNuclearCrossSection.hh"
#include "G4GammaNuclearXS.hh"

#include "G4HadronicParameters.hh"
#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"
#include "G4CrossSectionDataSetRegistry.hh"
 
// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmExtraPhysics);

//////////////////////////////////////

G4EmExtraPhysics::G4EmExtraPhysics(G4int ver): 
  G4VPhysicsConstructor("G4GammaLeptoNuclearPhys"),
  fGNLowEnergyLimit(200*CLHEP::MeV),
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
  // LEND cannot be used with low-energy model
  if(val) { fGNLowEnergyLimit = 0.0; }
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

void G4EmExtraPhysics::MuonToMuMu(G4bool val)
{
  mmumuActivated = val;
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

void G4EmExtraPhysics::SetUseGammaNuclearXS(G4bool val)
{
  fUseGammaNuclearXS = val;
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

void G4EmExtraPhysics::GammaNuclearLEModelLimit(G4double val)
{
  if(val <= CLHEP::MeV) { 
    fGNLowEnergyLimit = 0.0;

    // lowenergy model should not be applied at high energy
  } else if(val <= CLHEP::GeV) {
    fGNLowEnergyLimit = val;
    gLENDActivated = false;
  }
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

  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4LossTableManager* emManager = G4LossTableManager::Instance();

  if(gnActivated) { ConstructGammaElectroNuclear(); }

  if(munActivated) {
    G4MuonNuclearProcess* muNucProcess = new G4MuonNuclearProcess();
    G4MuonVDNuclearModel* muNucModel = new G4MuonVDNuclearModel();
    muNucProcess->RegisterMe(muNucModel);
    ph->RegisterProcess( muNucProcess, muonplus);
    ph->RegisterProcess( muNucProcess, muonminus);
  }
  if(gmumuActivated) {
    G4GammaConversionToMuons* theGammaToMuMu = new G4GammaConversionToMuons();
    theGammaToMuMu->SetCrossSecFactor(gmumuFactor);
    G4GammaGeneralProcess* sp = 
      static_cast<G4GammaGeneralProcess*>(emManager->GetGammaGeneralProcess());
    if(nullptr != sp) {
      sp->AddMMProcess(theGammaToMuMu);
    } else {
      ph->RegisterProcess(theGammaToMuMu, gamma);
    }
  }  
  if(mmumuActivated) {
    auto proc = new G4MuonToMuonPairProduction();
    ph->RegisterProcess(proc, muonplus);
    ph->RegisterProcess(proc, muonminus);
  }
  if(pmumuActivated) {
    G4AnnihiToMuPair* thePosiToMuMu = new G4AnnihiToMuPair();
    thePosiToMuMu->SetCrossSecFactor(pmumuFactor);
    ph->RegisterProcess(thePosiToMuMu, positron);
    G4AnnihiToMuPair* thePosiToTauTau = new G4AnnihiToMuPair("AnnihiToTauPair");
    thePosiToTauTau->SetCrossSecFactor(pmumuFactor);
    ph->RegisterProcess(thePosiToTauTau, positron);
  }  
  if(phadActivated) {
    G4eeToHadrons* thePosiToHadrons = new G4eeToHadrons();
    thePosiToHadrons->SetCrossSecFactor(phadFactor);
    ph->RegisterProcess(thePosiToHadrons, positron);
  }  
  if(synActivated) {
    G4SynchrotronRadiation* theSynchRad = new G4SynchrotronRadiation();
    ph->RegisterProcess( theSynchRad, electron);
    ph->RegisterProcess( theSynchRad, positron);
    if(synActivatedForAll) {
      ph->RegisterProcess( theSynchRad, muonplus);
      ph->RegisterProcess( theSynchRad, muonminus);

      ph->RegisterProcess( theSynchRad, G4Proton::Proton());
      ph->RegisterProcess( theSynchRad, G4AntiProton::AntiProton());
      ph->RegisterProcess( theSynchRad, G4PionPlus::PionPlus());
      ph->RegisterProcess( theSynchRad, G4PionMinus::PionMinus());
      ph->RegisterProcess( theSynchRad, G4GenericIon::GenericIon());
    }
  }
  if( fNuActivated )
  {
    G4ParticleDefinition* anuelectron = G4AntiNeutrinoE::AntiNeutrinoE();
    G4ParticleDefinition* nuelectron  = G4NeutrinoE::NeutrinoE(); 
    G4ParticleDefinition* anumuon     = G4AntiNeutrinoMu::AntiNeutrinoMu();
    G4ParticleDefinition* numuon      = G4NeutrinoMu::NeutrinoMu();
    G4ParticleDefinition* anutau      = G4AntiNeutrinoTau::AntiNeutrinoTau();
    G4ParticleDefinition* nutau       = G4NeutrinoTau::NeutrinoTau();

    G4NeutrinoElectronProcess* theNuEleProcess = 
      new G4NeutrinoElectronProcess(fNuDetectorName);
    G4NeutrinoElectronTotXsc* theNuEleTotXsc = new G4NeutrinoElectronTotXsc();

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

    // nu_mu nucleus interactions

    G4MuNeutrinoNucleusProcess* theNuMuNucleusProcess = new G4MuNeutrinoNucleusProcess(fNuDetectorName);
    G4MuNeutrinoNucleusTotXsc* theNuMuNucleusTotXsc = new G4MuNeutrinoNucleusTotXsc();
    
    if(fNuETotXscActivated)
    {
      theNuMuNucleusProcess->SetBiasingFactor(fNuNucleusBias);
    }
    theNuMuNucleusProcess->AddDataSet(theNuMuNucleusTotXsc);

    G4NuMuNucleusCcModel* numunuclcc = new G4NuMuNucleusCcModel();
    G4NuMuNucleusNcModel* numunuclnc = new G4NuMuNucleusNcModel();
    G4ANuMuNucleusCcModel* anumunuclcc = new G4ANuMuNucleusCcModel();
    G4ANuMuNucleusNcModel* anumunuclnc = new G4ANuMuNucleusNcModel();
    
    theNuMuNucleusProcess->RegisterMe(numunuclcc);
    theNuMuNucleusProcess->RegisterMe(numunuclnc);
    theNuMuNucleusProcess->RegisterMe(anumunuclcc);
    theNuMuNucleusProcess->RegisterMe(anumunuclnc);

    ph->RegisterProcess(theNuMuNucleusProcess, anumuon);
    ph->RegisterProcess(theNuMuNucleusProcess, numuon);    

    // nu_tau nucleus interactions

    G4TauNeutrinoNucleusProcess* theNuTauNucleusProcess = new G4TauNeutrinoNucleusProcess(fNuDetectorName);
    G4TauNeutrinoNucleusTotXsc* theNuTauNucleusTotXsc = new G4TauNeutrinoNucleusTotXsc();
    
    if(fNuETotXscActivated)
    {
      theNuTauNucleusProcess->SetBiasingFactor(fNuNucleusBias);
    }
    theNuTauNucleusProcess->AddDataSet(theNuTauNucleusTotXsc);

    G4NuTauNucleusCcModel* nutaunuclcc = new G4NuTauNucleusCcModel();
    G4NuTauNucleusNcModel* nutaunuclnc = new G4NuTauNucleusNcModel();
    G4ANuTauNucleusCcModel* anutaunuclcc = new G4ANuTauNucleusCcModel();
    G4ANuTauNucleusNcModel* anutaunuclnc = new G4ANuTauNucleusNcModel();
    
    theNuTauNucleusProcess->RegisterMe(nutaunuclcc);
    theNuTauNucleusProcess->RegisterMe(nutaunuclnc);
    theNuTauNucleusProcess->RegisterMe(anutaunuclcc);
    theNuTauNucleusProcess->RegisterMe(anutaunuclnc);

    ph->RegisterProcess(theNuTauNucleusProcess, anutau);
    ph->RegisterProcess(theNuTauNucleusProcess, nutau);    

    // nu_e nucleus interactions

    G4ElNeutrinoNucleusProcess* theNuElNucleusProcess = new G4ElNeutrinoNucleusProcess(fNuDetectorName);
    G4ElNeutrinoNucleusTotXsc* theNuElNucleusTotXsc = new G4ElNeutrinoNucleusTotXsc();
    
    if(fNuETotXscActivated)
    {
      theNuElNucleusProcess->SetBiasingFactor(fNuNucleusBias);
    }
    theNuElNucleusProcess->AddDataSet(theNuElNucleusTotXsc);

    G4NuElNucleusCcModel* nuelnuclcc = new G4NuElNucleusCcModel();
    G4NuElNucleusNcModel* nuelnuclnc = new G4NuElNucleusNcModel();
    G4ANuElNucleusCcModel* anuelnuclcc = new G4ANuElNucleusCcModel();
    G4ANuElNucleusNcModel* anuelnuclnc = new G4ANuElNucleusNcModel();
    
    theNuElNucleusProcess->RegisterMe(nuelnuclcc);
    theNuElNucleusProcess->RegisterMe(nuelnuclnc);
    theNuElNucleusProcess->RegisterMe(anuelnuclcc);
    theNuElNucleusProcess->RegisterMe(anuelnuclnc);

    ph->RegisterProcess(theNuElNucleusProcess, anuelectron);
    ph->RegisterProcess(theNuElNucleusProcess, nuelectron);    
  }
}

void G4EmExtraPhysics::ConstructGammaElectroNuclear()
{
  G4LossTableManager* emManager  = G4LossTableManager::Instance();
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  G4HadronInelasticProcess* gnuc = new G4HadronInelasticProcess( "photonNuclear", G4Gamma::Gamma() );
  auto xsreg = G4CrossSectionDataSetRegistry::Instance();
  G4VCrossSectionDataSet* xs = nullptr;
  if(fUseGammaNuclearXS) {
    xs = xsreg->GetCrossSectionDataSet("GammaNuclearXS");
    if(nullptr == xs) xs = new G4GammaNuclearXS();
  } else {
    xs = xsreg->GetCrossSectionDataSet("PhotoNuclearXS");
    if(nullptr == xs) xs = new G4PhotoNuclearCrossSection(); 
  }
  gnuc->AddDataSet(xs);

  G4QGSModel< G4GammaParticipants >* theStringModel = 
    new G4QGSModel< G4GammaParticipants >;
  G4QGSMFragmentation* theFrag = new G4QGSMFragmentation();
  G4ExcitedStringDecay* theStringDecay = new G4ExcitedStringDecay(theFrag);
  theStringModel->SetFragmentationModel(theStringDecay);

  G4GeneratorPrecompoundInterface* theCascade = 
    new G4GeneratorPrecompoundInterface;

  G4TheoFSGenerator* theModel = new G4TheoFSGenerator();
  theModel->SetTransport(theCascade);
  theModel->SetHighEnergyGenerator(theStringModel);

  G4HadronicParameters* param = G4HadronicParameters::Instance();

  G4CascadeInterface* cascade = new G4CascadeInterface;

  // added low-energy model LEND disabled
  if (fGNLowEnergyLimit > 0.0) { 
    G4LowEGammaNuclearModel* lemod = new G4LowEGammaNuclearModel();
    lemod->SetMaxEnergy(fGNLowEnergyLimit);
    gnuc->RegisterMe(lemod);
    cascade->SetMinEnergy(fGNLowEnergyLimit - CLHEP::MeV);
  }
  cascade->SetMaxEnergy(param->GetMaxEnergyTransitionFTF_Cascade());
  gnuc->RegisterMe(cascade);
  theModel->SetMinEnergy(param->GetMinEnergyTransitionFTF_Cascade());
  theModel->SetMaxEnergy(param->GetMaxEnergy());
  gnuc->RegisterMe(theModel);
  
  G4GammaGeneralProcess* gproc = 
    (G4GammaGeneralProcess*)emManager->GetGammaGeneralProcess();
  if(gproc != nullptr) {
    gproc->AddHadProcess(gnuc);
  } else {
    // LEND may be activated if the general process is not activated
    ph->RegisterProcess(gnuc, G4Gamma::Gamma());
    if(gLENDActivated) { ConstructLENDGammaNuclear(cascade, gnuc); }
  }

  if(eActivated) {
    G4ElectronNuclearProcess* enuc = new G4ElectronNuclearProcess;
    G4PositronNuclearProcess* pnuc = new G4PositronNuclearProcess;
    G4ElectroVDNuclearModel* eModel = new G4ElectroVDNuclearModel;

    enuc->RegisterMe(eModel);
    pnuc->RegisterMe(eModel);

    G4GammaGeneralProcess* eproc = 
      (G4GammaGeneralProcess*)emManager->GetElectronGeneralProcess();
    if(eproc != nullptr) {
      eproc->AddHadProcess(enuc);
    } else {
      ph->RegisterProcess(enuc, G4Electron::Electron());
    }
  
    G4GammaGeneralProcess* pproc = 
      (G4GammaGeneralProcess*)emManager->GetPositronGeneralProcess();
    if(pproc != nullptr) {
      pproc->AddHadProcess(pnuc);
    } else {
      ph->RegisterProcess(pnuc, G4Positron::Positron());
    }
  }
}

void G4EmExtraPhysics::ConstructLENDGammaNuclear(
     G4CascadeInterface* cascade, G4HadronInelasticProcess* gnuc)
{
  if (G4FindDataDir("G4LENDDATA") == nullptr ) {
    G4String message = "\n Skipping activation of Low Energy Nuclear Data (LEND) model for gamma nuclear interactions.\n The LEND model needs data files and they are available from ftp://gdo-nuclear.ucllnl.org/GND_after2013/GND_v1.3.tar.gz.\n Please set the environment variable G4LENDDATA to point to the directory named v1.3 extracted from the archive file.\n"; 
    G4Exception( "G4EmExtraPhysics::ConstructLENDGammaNuclear()"
                 , "G4LENDBertiniGammaElectroNuclearBuilder001"
                 , JustWarning , message);
    return;
  }
   
  cascade->SetMinEnergy(19.9*MeV);
  G4LENDorBERTModel* theGammaReactionLowE = 
    new G4LENDorBERTModel( G4Gamma::Gamma() );
  theGammaReactionLowE->DumpLENDTargetInfo(true);
  G4LENDCombinedCrossSection* theGammaCrossSectionLowE = 
    new G4LENDCombinedCrossSection( G4Gamma::Gamma() );
  theGammaReactionLowE->SetMaxEnergy(20*MeV);
  gnuc->RegisterMe(theGammaReactionLowE);
  gnuc->AddDataSet(theGammaCrossSectionLowE);   
}
