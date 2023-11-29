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
/// \file HadronicGenerator.cc
/// \brief Implementation of the HadronicGenerator class
//
//------------------------------------------------------------------------
// Class: HadronicGenerator
// Author: Alberto Ribon (CERN EP/SFT)
// Date: May 2020
//------------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HadronicGenerator.hh"
#include <iomanip>
#include "globals.hh"
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4Material.hh"
#include "G4ProcessManager.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayPhysics.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4StateManager.hh"
#include "G4TouchableHistory.hh"
#include "G4TransportationManager.hh"

#include "G4PionMinus.hh"
#include "G4PionPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonZeroLong.hh"
#include "G4KaonZeroShort.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4Lambda.hh"
#include "G4SigmaPlus.hh"
#include "G4SigmaZero.hh"
#include "G4SigmaMinus.hh"
#include "G4XiMinus.hh"
#include "G4XiZero.hh"
#include "G4OmegaMinus.hh"
#include "G4AntiProton.hh"
#include "G4AntiNeutron.hh"
#include "G4AntiDeuteron.hh"
#include "G4AntiTriton.hh"
#include "G4AntiHe3.hh"
#include "G4AntiAlpha.hh"
#include "G4AntiLambda.hh"
#include "G4AntiSigmaPlus.hh"
#include "G4AntiSigmaZero.hh"
#include "G4AntiSigmaMinus.hh"
#include "G4AntiXiMinus.hh"
#include "G4AntiXiZero.hh"
#include "G4AntiOmegaMinus.hh"
#include "G4GenericIon.hh"

#include "G4DMesonPlus.hh"
#include "G4DMesonMinus.hh"
#include "G4DMesonZero.hh"
#include "G4AntiDMesonZero.hh"
#include "G4DsMesonPlus.hh"
#include "G4DsMesonMinus.hh"
#include "G4BMesonPlus.hh"
#include "G4BMesonMinus.hh"
#include "G4BMesonZero.hh"
#include "G4AntiBMesonZero.hh"
#include "G4BsMesonZero.hh"
#include "G4AntiBsMesonZero.hh"
#include "G4BcMesonPlus.hh"
#include "G4BcMesonMinus.hh"
#include "G4LambdacPlus.hh"
#include "G4AntiLambdacPlus.hh"
#include "G4XicPlus.hh"
#include "G4AntiXicPlus.hh"
#include "G4XicZero.hh"
#include "G4AntiXicZero.hh"
#include "G4OmegacZero.hh"
#include "G4AntiOmegacZero.hh"
#include "G4Lambdab.hh"
#include "G4AntiLambdab.hh"
#include "G4XibZero.hh"
#include "G4AntiXibZero.hh"
#include "G4XibMinus.hh"
#include "G4AntiXibMinus.hh"
#include "G4OmegabMinus.hh"
#include "G4AntiOmegabMinus.hh"

#include "G4HyperTriton.hh"
#include "G4AntiHyperTriton.hh"
#include "G4HyperAlpha.hh"
#include "G4AntiHyperAlpha.hh"
#include "G4HyperH4.hh"
#include "G4AntiHyperH4.hh"
#include "G4DoubleHyperH4.hh"
#include "G4AntiDoubleHyperH4.hh"
#include "G4DoubleHyperDoubleNeutron.hh"
#include "G4AntiDoubleHyperDoubleNeutron.hh"
#include "G4HyperHe5.hh"
#include "G4AntiHyperHe5.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4HadronicParameters.hh"

#include "G4CascadeInterface.hh"
#include "G4TheoFSGenerator.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4FTFModel.hh"
#include "G4BinaryCascade.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4INCLXXInterface.hh"
#include "G4AblaInterface.hh"
#include "G4QuasiElasticChannel.hh"
#include "G4QGSMFragmentation.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"

#include "G4VCrossSectionDataSet.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4BGGPionInelasticXS.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4ChipsHyperonInelasticXS.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4ComponentGGNuclNuclXsc.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HadronicGenerator::HadronicGenerator( const G4String physicsCase ) :
  fPhysicsCase( physicsCase ), fPhysicsCaseIsSupported( false ),
  fLastHadronicProcess( nullptr ), fPartTable( nullptr )
{
  // The constructor set-ups all the particles, models, cross sections and
  // hadronic inelastic processes.
  // This should be done only once for each application.
  // In the case of a multi-threaded application using this class,
  // the constructor should be invoked for each thread,
  // i.e. one instance of the class should be kept per thread.
  // The particles and processes that are created in this constructor
  // will then be used by the method GenerateInteraction at each interaction.
  // Notes:
  // - Neither the hadronic models nor the cross sections are used directly
  //   by the method GenerateInteraction, but they are associated to the
  //   hadronic processes and used by Geant4 to simulate the collision;
  // - Although the class generates only final states, but not free mean paths,
  //   inelastic hadron-nuclear cross sections are needed by Geant4 to sample
  //   the target nucleus from the target material.
  
  // Definition of particles
  G4GenericIon* gion = G4GenericIon::Definition();
  gion->SetProcessManager( new G4ProcessManager( gion ) );
  G4DecayPhysics* decays = new G4DecayPhysics;
  decays->ConstructParticle();  
  fPartTable = G4ParticleTable::GetParticleTable();
  fPartTable->SetReadiness();
  G4IonTable* ions = fPartTable->GetIonTable();
  ions->CreateAllIon();
  ions->CreateAllIsomer();

  // Build BERT model
  G4CascadeInterface* theBERTmodel = new G4CascadeInterface;

  // Build BIC model
  G4BinaryCascade* theBICmodel = new G4BinaryCascade;
  G4PreCompoundModel* thePreEquilib = new G4PreCompoundModel( new G4ExcitationHandler );
  theBICmodel->SetDeExcitation( thePreEquilib );

  // Build BinaryLightIon model
  G4PreCompoundModel* thePreEquilibBis = new G4PreCompoundModel( new G4ExcitationHandler );
  G4BinaryLightIonReaction* theIonBICmodel = new G4BinaryLightIonReaction( thePreEquilibBis );

  // Build the INCL model
  G4INCLXXInterface* theINCLmodel = new G4INCLXXInterface;
  const G4bool useAblaDeExcitation = false;  // By default INCL uses Preco: set "true" to use
                                             // ABLA DeExcitation
  if ( theINCLmodel && useAblaDeExcitation ) {
    G4AblaInterface* theAblaInterface = new G4AblaInterface;
    theINCLmodel->SetDeExcitation( theAblaInterface );
  }
  
  // Build the FTFP model (FTF/Preco) : 4 instances with different kinetic energy intervals.
  // (Notice that these kinetic energy intervals are applied per nucleons, so they are fine
  // for all types of hadron and ion projectile).
  // Model instance without energy constraint.
  // (Used for the case of FTFP model, and for light anti-ions in all physics lists.)
  G4TheoFSGenerator* theFTFPmodel = new G4TheoFSGenerator( "FTFP" );
  theFTFPmodel->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  G4GeneratorPrecompoundInterface* theCascade = new G4GeneratorPrecompoundInterface;
  theCascade->SetDeExcitation( thePreEquilib );
  theFTFPmodel->SetTransport( theCascade );
  G4LundStringFragmentation* theLundFragmentation = new G4LundStringFragmentation;
  G4ExcitedStringDecay* theStringDecay = new G4ExcitedStringDecay( theLundFragmentation );
  G4FTFModel* theStringModel = new G4FTFModel;
  theStringModel->SetFragmentationModel( theStringDecay );

  // If the following line is set, then the square of the impact parameter is sampled
  // randomly from a flat distribution in the range [ Bmin*Bmin, Bmax*Bmax ]
  //***LOOKHERE*** CHOOSE IMPACT PARAMETER MIN & MAX
  //theStringModel->SetBminBmax( 0.0, 2.0*fermi );

  theFTFPmodel->SetHighEnergyGenerator( theStringModel );
  // Model instance with constraint to be above a kinetic energy threshold.
  // (Used for ions in all physics lists, and, in the case of non-QGS-based physics lists,
  // also for pions, kaons, nucleons and hyperons.)
  G4TheoFSGenerator* theFTFPmodel_aboveThreshold = new G4TheoFSGenerator( "FTFP" );
  theFTFPmodel_aboveThreshold->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  theFTFPmodel_aboveThreshold->SetTransport( theCascade );
  theFTFPmodel_aboveThreshold->SetHighEnergyGenerator( theStringModel );
  // Model instance with constraint to be within two kinetic energy thresholds.
  // (Used in the case of QGS-based physics lists for pions, kaons, nucleons and hyperons.)
  G4TheoFSGenerator* theFTFPmodel_constrained = new G4TheoFSGenerator( "FTFP" );
  theFTFPmodel_constrained->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  theFTFPmodel_constrained->SetTransport( theCascade );
  theFTFPmodel_constrained->SetHighEnergyGenerator( theStringModel );
  // Model instance to be used down to zero kinetic energy, with eventual constraint
  // - in the case of QGS-based physics lists - to be below a kinetic energy threshold.
  // (Used for anti-baryons, anti-hyperons, and charmed and bottom hadrons.)
  G4TheoFSGenerator* theFTFPmodel_belowThreshold = new G4TheoFSGenerator( "FTFP" );
  theFTFPmodel_belowThreshold->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  theFTFPmodel_belowThreshold->SetTransport( theCascade );
  theFTFPmodel_belowThreshold->SetHighEnergyGenerator( theStringModel );

  // Build the QGSP model (QGS/Preco)
  G4TheoFSGenerator* theQGSPmodel = new G4TheoFSGenerator( "QGSP" );
  theQGSPmodel->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  theQGSPmodel->SetTransport( theCascade );
  G4QGSMFragmentation* theQgsmFragmentation = new G4QGSMFragmentation;
  G4ExcitedStringDecay* theQgsmStringDecay = new G4ExcitedStringDecay( theQgsmFragmentation );
  G4VPartonStringModel* theQgsmStringModel = new G4QGSModel< G4QGSParticipants >;
  theQgsmStringModel->SetFragmentationModel( theQgsmStringDecay );
  theQGSPmodel->SetHighEnergyGenerator( theQgsmStringModel );
  G4QuasiElasticChannel* theQuasiElastic = new G4QuasiElasticChannel;  // QGSP uses quasi-elastic
  theQGSPmodel->SetQuasiElasticChannel( theQuasiElastic );

  // For the case of "physics-list proxies", select the energy range for each hadronic model.
  // Note: the transition energy between hadronic models vary between physics lists,
  //       type of hadrons, and version of Geant4. Here, for simplicity, we use an uniform
  //       energy transition for all types of hadrons and regardless of the Geant4 version;
  //       moreover, for "FTFP_INCLXX" we use a different energy transition range
  //       between FTFP and INCL than in the real physics list.
  if ( fPhysicsCase == "FTFP_BERT_ATL"  ||
       fPhysicsCase == "FTFP_BERT"      ||
       fPhysicsCase == "FTFP_INCLXX"    ||
       fPhysicsCase == "QGSP_BERT"      ||
       fPhysicsCase == "QGSP_BIC" ) {
    const G4double ftfpMinE = G4HadronicParameters::Instance()->GetMinEnergyTransitionFTF_Cascade();
    const G4double bertMaxE = G4HadronicParameters::Instance()->GetMaxEnergyTransitionFTF_Cascade();
    const G4double ftfpMinE_ATL =  9.0*CLHEP::GeV;
    const G4double bertMaxE_ATL = 12.0*CLHEP::GeV;
    const G4double ftfpMaxE = G4HadronicParameters::Instance()->GetMaxEnergyTransitionQGS_FTF();
    const G4double qgspMinE = G4HadronicParameters::Instance()->GetMinEnergyTransitionQGS_FTF();
    theFTFPmodel->SetMinEnergy( 0.0 );
    theFTFPmodel_belowThreshold->SetMinEnergy( 0.0 );
    if ( fPhysicsCase == "FTFP_BERT_ATL" ) {
      theBERTmodel->SetMaxEnergy( bertMaxE_ATL );
      theIonBICmodel->SetMaxEnergy( bertMaxE_ATL );
      theFTFPmodel_aboveThreshold->SetMinEnergy( ftfpMinE_ATL );
      theFTFPmodel_constrained->SetMinEnergy( ftfpMinE_ATL );
    } else {
      theBERTmodel->SetMaxEnergy( bertMaxE );
      theIonBICmodel->SetMaxEnergy( bertMaxE );
      theFTFPmodel_aboveThreshold->SetMinEnergy( ftfpMinE );
      theFTFPmodel_constrained->SetMinEnergy( ftfpMinE );
    }
    if ( fPhysicsCase == "FTFP_INCLXX" ) {
      theINCLmodel->SetMaxEnergy( bertMaxE );
    }
    if ( fPhysicsCase == "QGSP_BERT"  ||
         fPhysicsCase == "QGSP_BIC" ) {
      theFTFPmodel_constrained->SetMaxEnergy( ftfpMaxE );
      theFTFPmodel_belowThreshold->SetMaxEnergy( ftfpMaxE );
      theQGSPmodel->SetMinEnergy( qgspMinE );
      theBICmodel->SetMaxEnergy( bertMaxE );
    }
  }

  // Cross sections (needed by Geant4 to sample the target nucleus from the target material)
  G4VCrossSectionDataSet* thePionMinusXSdata =
    new G4BGGPionInelasticXS( G4PionMinus::Definition() );
  thePionMinusXSdata->BuildPhysicsTable( *(G4PionMinus::Definition()) );
  G4VCrossSectionDataSet* thePionPlusXSdata =
    new G4BGGPionInelasticXS( G4PionPlus::Definition() );
  thePionPlusXSdata->BuildPhysicsTable( *(G4PionPlus::Definition()) );
  G4VCrossSectionDataSet* theKaonXSdata =
    new G4CrossSectionInelastic( new G4ComponentGGHadronNucleusXsc );
  theKaonXSdata->BuildPhysicsTable( *(G4KaonMinus::Definition()) );
  theKaonXSdata->BuildPhysicsTable( *(G4KaonPlus::Definition()) );
  theKaonXSdata->BuildPhysicsTable( *(G4KaonZeroLong::Definition()) );
  theKaonXSdata->BuildPhysicsTable( *(G4KaonZeroShort::Definition()) );
  G4VCrossSectionDataSet* theProtonXSdata = new G4BGGNucleonInelasticXS( G4Proton::Proton() );
  theProtonXSdata->BuildPhysicsTable( *(G4Proton::Definition()) );
  G4VCrossSectionDataSet* theNeutronXSdata = new G4NeutronInelasticXS;
  theNeutronXSdata->BuildPhysicsTable( *(G4Neutron::Definition()) );
  // For hyperon and anti-hyperons we can use either Chips or, for G4 >= 10.5,
  // Glauber-Gribov cross sections
  //G4VCrossSectionDataSet* theHyperonsXSdata = new G4ChipsHyperonInelasticXS;
  G4VCrossSectionDataSet* theHyperonsXSdata =
    new G4CrossSectionInelastic( new G4ComponentGGHadronNucleusXsc );
  G4VCrossSectionDataSet* theAntibaryonsXSdata =
    new G4CrossSectionInelastic( new G4ComponentAntiNuclNuclearXS );
  G4VCrossSectionDataSet* theNuclNuclXSdata =
    new G4CrossSectionInelastic( new G4ComponentGGNuclNuclXsc );

  // Set up inelastic processes : store them in a map (with particle definition as key)
  //                              for convenience
  typedef std::pair< G4ParticleDefinition*, G4HadronicProcess* > ProcessPair;
  G4HadronicProcess* thePionMinusInelasticProcess =
    new G4HadronInelasticProcess( "pi-Inelastic", G4PionMinus::Definition() );    
  fProcessMap.insert( ProcessPair( G4PionMinus::Definition(), thePionMinusInelasticProcess ) );
  G4HadronicProcess* thePionPlusInelasticProcess =
    new G4HadronInelasticProcess( "pi+Inelastic", G4PionPlus::Definition() );
  fProcessMap.insert( ProcessPair( G4PionPlus::Definition(), thePionPlusInelasticProcess ) );
  G4HadronicProcess* theKaonMinusInelasticProcess =
    new G4HadronInelasticProcess( "kaon-Inelastic", G4KaonMinus::Definition() );
  fProcessMap.insert( ProcessPair( G4KaonMinus::Definition(), theKaonMinusInelasticProcess ) );
  G4HadronicProcess* theKaonPlusInelasticProcess =
    new G4HadronInelasticProcess( "kaon+Inelastic", G4KaonPlus::Definition() );
  fProcessMap.insert( ProcessPair( G4KaonPlus::Definition(), theKaonPlusInelasticProcess ) );
  G4HadronicProcess* theKaonZeroLInelasticProcess =
    new G4HadronInelasticProcess( "kaon0LInelastic", G4KaonZeroLong::Definition() );
  fProcessMap.insert( ProcessPair( G4KaonZeroLong::Definition(), theKaonZeroLInelasticProcess ) );
  G4HadronicProcess* theKaonZeroSInelasticProcess =
    new G4HadronInelasticProcess( "kaon0SInelastic", G4KaonZeroShort::Definition() );
  fProcessMap.insert( ProcessPair( G4KaonZeroShort::Definition(), theKaonZeroSInelasticProcess ) );
  G4HadronicProcess* theProtonInelasticProcess =
    new G4HadronInelasticProcess( "protonInelastic", G4Proton::Definition() );
  fProcessMap.insert( ProcessPair( G4Proton::Definition(), theProtonInelasticProcess ) );
  G4HadronicProcess* theNeutronInelasticProcess =
    new G4HadronInelasticProcess( "neutronInelastic", G4Neutron::Definition() );
  fProcessMap.insert( ProcessPair( G4Neutron::Definition(), theNeutronInelasticProcess ) );
  G4HadronicProcess* theDeuteronInelasticProcess =
    new G4HadronInelasticProcess( "dInelastic", G4Deuteron::Definition() );
  fProcessMap.insert( ProcessPair( G4Deuteron::Definition(), theDeuteronInelasticProcess ) );
  G4HadronicProcess* theTritonInelasticProcess =
    new G4HadronInelasticProcess( "tInelastic", G4Triton::Definition() );
  fProcessMap.insert( ProcessPair( G4Triton::Definition(), theTritonInelasticProcess ) );
  G4HadronicProcess* theHe3InelasticProcess =
    new G4HadronInelasticProcess( "he3Inelastic", G4He3::Definition() );
  fProcessMap.insert( ProcessPair( G4He3::Definition(), theHe3InelasticProcess ) );
  G4HadronicProcess* theAlphaInelasticProcess =
    new G4HadronInelasticProcess( "alphaInelastic", G4Alpha::Definition() ); 
  fProcessMap.insert( ProcessPair( G4Alpha::Definition(), theAlphaInelasticProcess ) );
  G4HadronicProcess* theIonInelasticProcess =
    new G4HadronInelasticProcess( "ionInelastic", G4GenericIon::Definition() );
  fProcessMap.insert( ProcessPair( G4GenericIon::Definition(), theIonInelasticProcess ) );
  G4HadronicProcess* theLambdaInelasticProcess =
    new G4HadronInelasticProcess( "lambdaInelastic", G4Lambda::Definition() );
  fProcessMap.insert( ProcessPair( G4Lambda::Definition(), theLambdaInelasticProcess ) );
  G4HadronicProcess* theSigmaMinusInelasticProcess =
    new G4HadronInelasticProcess( "sigma-Inelastic", G4SigmaMinus::Definition() );
  fProcessMap.insert( ProcessPair( G4SigmaMinus::Definition(), theSigmaMinusInelasticProcess ) );
  G4HadronicProcess* theSigmaPlusInelasticProcess =
    new G4HadronInelasticProcess( "sigma+Inelastic", G4SigmaPlus::Definition() );
  fProcessMap.insert( ProcessPair( G4SigmaPlus::Definition(), theSigmaPlusInelasticProcess ) );
  G4HadronicProcess* theXiMinusInelasticProcess =
    new G4HadronInelasticProcess( "xi-Inelastic", G4XiMinus::Definition() );
  fProcessMap.insert( ProcessPair( G4XiMinus::Definition(), theXiMinusInelasticProcess ) );
  G4HadronicProcess* theXiZeroInelasticProcess =
    new G4HadronInelasticProcess( "xi0Inelastic", G4XiZero::Definition() );
  fProcessMap.insert( ProcessPair( G4XiZero::Definition(), theXiZeroInelasticProcess ) );
  G4HadronicProcess* theOmegaMinusInelasticProcess =
    new G4HadronInelasticProcess( "omega-Inelastic", G4OmegaMinus::Definition() );
  fProcessMap.insert( ProcessPair( G4OmegaMinus::Definition(), theOmegaMinusInelasticProcess ) );
  G4HadronicProcess* theAntiProtonInelasticProcess =
    new G4HadronInelasticProcess( "anti_protonInelastic", G4AntiProton::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiProton::Definition(), theAntiProtonInelasticProcess ) );
  G4HadronicProcess* theAntiNeutronInelasticProcess =
    new G4HadronInelasticProcess( "anti_neutronInelastic", G4AntiNeutron::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiNeutron::Definition(), theAntiNeutronInelasticProcess ) );
  G4HadronicProcess* theAntiDeuteronInelasticProcess =
    new G4HadronInelasticProcess( "anti_deuteronInelastic", G4AntiDeuteron::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiDeuteron::Definition(),
                                   theAntiDeuteronInelasticProcess ) );
  G4HadronicProcess* theAntiTritonInelasticProcess =
    new G4HadronInelasticProcess( "anti_tritonInelastic", G4AntiTriton::Definition() );    
  fProcessMap.insert( ProcessPair( G4AntiTriton::Definition(), theAntiTritonInelasticProcess ) );
  G4HadronicProcess* theAntiHe3InelasticProcess =
    new G4HadronInelasticProcess( "anti_He3Inelastic", G4AntiHe3::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiHe3::Definition(), theAntiHe3InelasticProcess ) );
  G4HadronicProcess* theAntiAlphaInelasticProcess =
    new G4HadronInelasticProcess( "anti_alphaInelastic", G4AntiAlpha::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiAlpha::Definition(), theAntiAlphaInelasticProcess ) );
  G4HadronicProcess* theAntiLambdaInelasticProcess =
    new G4HadronInelasticProcess( "anti-lambdaInelastic", G4AntiLambda::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiLambda::Definition(), theAntiLambdaInelasticProcess ) );
  G4HadronicProcess* theAntiSigmaMinusInelasticProcess =
    new G4HadronInelasticProcess( "anti_sigma-Inelastic", G4AntiSigmaMinus::Definition() );  
  fProcessMap.insert( ProcessPair( G4AntiSigmaMinus::Definition(),
                                   theAntiSigmaMinusInelasticProcess ) );
  G4HadronicProcess* theAntiSigmaPlusInelasticProcess =
    new G4HadronInelasticProcess( "anti_sigma+Inelastic", G4AntiSigmaPlus::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiSigmaPlus::Definition(),
                                   theAntiSigmaPlusInelasticProcess ) );
  G4HadronicProcess* theAntiXiMinusInelasticProcess =
    new G4HadronInelasticProcess( "anti_xi-Inelastic", G4AntiXiMinus::Definition() );   
  fProcessMap.insert( ProcessPair( G4AntiXiMinus::Definition(), theAntiXiMinusInelasticProcess ) );
  G4HadronicProcess* theAntiXiZeroInelasticProcess =
    new G4HadronInelasticProcess( "anti_xi0Inelastic", G4AntiXiZero::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiXiZero::Definition(), theAntiXiZeroInelasticProcess ) );
  G4HadronicProcess* theAntiOmegaMinusInelasticProcess =
    new G4HadronInelasticProcess( "anti_omega-Inelastic", G4AntiOmegaMinus::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiOmegaMinus::Definition(),
                                   theAntiOmegaMinusInelasticProcess ) );

  G4HadronicProcess* theDPlusInelasticProcess =
    new G4HadronInelasticProcess( "D+Inelastic", G4DMesonPlus::Definition() );
  fProcessMap.insert( ProcessPair( G4DMesonPlus::Definition(), theDPlusInelasticProcess ) );
  G4HadronicProcess* theDMinusInelasticProcess =
    new G4HadronInelasticProcess( "D-Inelastic", G4DMesonMinus::Definition() );
  fProcessMap.insert( ProcessPair( G4DMesonMinus::Definition(), theDMinusInelasticProcess ) );
  G4HadronicProcess* theDZeroInelasticProcess =
    new G4HadronInelasticProcess( "D0Inelastic", G4DMesonZero::Definition() );
  fProcessMap.insert( ProcessPair( G4DMesonZero::Definition(), theDZeroInelasticProcess ) );
  G4HadronicProcess* theAntiDZeroInelasticProcess =
    new G4HadronInelasticProcess( "anti_D0Inelastic", G4AntiDMesonZero::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiDMesonZero::Definition(),
                                   theAntiDZeroInelasticProcess ) );
  G4HadronicProcess* theDsPlusInelasticProcess =
    new G4HadronInelasticProcess( "Ds+Inelastic", G4DsMesonPlus::Definition() );
  fProcessMap.insert( ProcessPair( G4DsMesonPlus::Definition(), theDsPlusInelasticProcess ) );
  G4HadronicProcess* theDsMinusInelasticProcess =
    new G4HadronInelasticProcess( "Ds-Inelastic", G4DsMesonMinus::Definition() );
  fProcessMap.insert( ProcessPair( G4DsMesonMinus::Definition(), theDsMinusInelasticProcess ) );
  G4HadronicProcess* theBPlusInelasticProcess =
    new G4HadronInelasticProcess( "B+Inelastic", G4BMesonPlus::Definition() );
  fProcessMap.insert( ProcessPair( G4BMesonPlus::Definition(), theBPlusInelasticProcess ) );
  G4HadronicProcess* theBMinusInelasticProcess =
    new G4HadronInelasticProcess( "B-Inelastic", G4BMesonMinus::Definition() );
  fProcessMap.insert( ProcessPair( G4BMesonMinus::Definition(), theBMinusInelasticProcess ) );
  G4HadronicProcess* theBZeroInelasticProcess =
    new G4HadronInelasticProcess( "B0Inelastic", G4BMesonZero::Definition() );
  fProcessMap.insert( ProcessPair( G4BMesonZero::Definition(), theBZeroInelasticProcess ) );
  G4HadronicProcess* theAntiBZeroInelasticProcess =
    new G4HadronInelasticProcess( "anti_B0Inelastic", G4AntiBMesonZero::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiBMesonZero::Definition(),
                                   theAntiBZeroInelasticProcess ) );
  G4HadronicProcess* theBsZeroInelasticProcess =
    new G4HadronInelasticProcess( "Bs0Inelastic", G4BsMesonZero::Definition() );
  fProcessMap.insert( ProcessPair( G4BsMesonZero::Definition(), theBsZeroInelasticProcess ) );
  G4HadronicProcess* theAntiBsZeroInelasticProcess =
    new G4HadronInelasticProcess( "anti_Bs0Inelastic", G4AntiBsMesonZero::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiBsMesonZero::Definition(),
                                   theAntiBsZeroInelasticProcess ) );
  G4HadronicProcess* theBcPlusInelasticProcess =
    new G4HadronInelasticProcess( "Bc+Inelastic", G4BcMesonPlus::Definition() );
  fProcessMap.insert( ProcessPair( G4BcMesonPlus::Definition(), theBcPlusInelasticProcess ) );
  G4HadronicProcess* theBcMinusInelasticProcess =
    new G4HadronInelasticProcess( "Bc-Inelastic", G4BcMesonMinus::Definition() );
  fProcessMap.insert( ProcessPair( G4BcMesonMinus::Definition(), theBcMinusInelasticProcess ) );
  G4HadronicProcess* theLambdacPlusInelasticProcess =
    new G4HadronInelasticProcess( "lambda_c+Inelastic", G4LambdacPlus::Definition() );
  fProcessMap.insert( ProcessPair( G4LambdacPlus::Definition(),
                                   theLambdacPlusInelasticProcess ) );
  G4HadronicProcess* theAntiLambdacPlusInelasticProcess =
    new G4HadronInelasticProcess( "anti_lambda_c+Inelastic", G4AntiLambdacPlus::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiLambdacPlus::Definition(),
                                   theAntiLambdacPlusInelasticProcess ) );
  G4HadronicProcess* theXicPlusInelasticProcess =
    new G4HadronInelasticProcess( "xi_c+Inelastic", G4XicPlus::Definition() );
  fProcessMap.insert( ProcessPair( G4XicPlus::Definition(), theXicPlusInelasticProcess ) );
  G4HadronicProcess* theAntiXicPlusInelasticProcess =
    new G4HadronInelasticProcess( "anti_xi_c+Inelastic", G4AntiXicPlus::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiXicPlus::Definition(),
                                   theAntiXicPlusInelasticProcess ) );
  G4HadronicProcess* theXicZeroInelasticProcess =
    new G4HadronInelasticProcess( "xi_c0Inelastic", G4XicZero::Definition() );
  fProcessMap.insert( ProcessPair( G4XicZero::Definition(), theXicZeroInelasticProcess ) );
  G4HadronicProcess* theAntiXicZeroInelasticProcess =
    new G4HadronInelasticProcess( "anti_xi_c0Inelastic", G4AntiXicZero::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiXicZero::Definition(),
                                   theAntiXicZeroInelasticProcess ) );
  G4HadronicProcess* theOmegacZeroInelasticProcess =
    new G4HadronInelasticProcess( "omega_c0Inelastic", G4OmegacZero::Definition() );
  fProcessMap.insert( ProcessPair( G4OmegacZero::Definition(), theOmegacZeroInelasticProcess ) );
  G4HadronicProcess* theAntiOmegacZeroInelasticProcess =
    new G4HadronInelasticProcess( "anti_omega_c0Inelastic", G4AntiOmegacZero::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiOmegacZero::Definition(),
                                   theAntiOmegacZeroInelasticProcess ) );
  G4HadronicProcess* theLambdabInelasticProcess =
    new G4HadronInelasticProcess( "lambda_bInelastic", G4Lambdab::Definition() );
  fProcessMap.insert( ProcessPair( G4Lambdab::Definition(), theLambdabInelasticProcess ) );
  G4HadronicProcess* theAntiLambdabInelasticProcess =
    new G4HadronInelasticProcess( "anti_lambda_bInelastic", G4AntiLambdab::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiLambdab::Definition(),
                                   theAntiLambdabInelasticProcess ) );
  G4HadronicProcess* theXibZeroInelasticProcess =
    new G4HadronInelasticProcess( "xi_b0Inelastic", G4XibZero::Definition() );
  fProcessMap.insert( ProcessPair( G4XibZero::Definition(), theXibZeroInelasticProcess ) );
  G4HadronicProcess* theAntiXibZeroInelasticProcess =
    new G4HadronInelasticProcess( "anti_xi_b0Inelastic", G4AntiXibZero::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiXibZero::Definition(),
                                   theAntiXibZeroInelasticProcess ) );
  G4HadronicProcess* theXibMinusInelasticProcess =
    new G4HadronInelasticProcess( "xi_b-Inelastic", G4XibMinus::Definition() );
  fProcessMap.insert( ProcessPair( G4XibMinus::Definition(), theXibMinusInelasticProcess ) );
  G4HadronicProcess* theAntiXibMinusInelasticProcess =
    new G4HadronInelasticProcess( "anti_xi_b-Inelastic", G4AntiXibMinus::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiXibMinus::Definition(),
                                   theAntiXibMinusInelasticProcess ) );
  G4HadronicProcess* theOmegabMinusInelasticProcess =
    new G4HadronInelasticProcess( "omega_b-Inelastic", G4OmegabMinus::Definition() );
  fProcessMap.insert( ProcessPair( G4OmegabMinus::Definition(),
                                   theOmegabMinusInelasticProcess ) );
  G4HadronicProcess* theAntiOmegabMinusInelasticProcess =
    new G4HadronInelasticProcess( "anti_omega_b-Inelastic", G4AntiOmegabMinus::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiOmegabMinus::Definition(),
                                   theAntiOmegabMinusInelasticProcess ) );

  G4HadronicProcess* theHyperTritonInelasticProcess =
    new G4HadronInelasticProcess( "hypertritonInelastic", G4HyperTriton::Definition() );
  fProcessMap.insert( ProcessPair( G4HyperTriton::Definition(),
                                   theHyperTritonInelasticProcess ) );
  G4HadronicProcess* theAntiHyperTritonInelasticProcess =
    new G4HadronInelasticProcess( "anti_hypertritonInelastic", G4AntiHyperTriton::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiHyperTriton::Definition(),
                                   theAntiHyperTritonInelasticProcess ) );
  G4HadronicProcess* theHyperAlphaInelasticProcess =
    new G4HadronInelasticProcess( "hyperalphaInelastic", G4HyperAlpha::Definition() );
  fProcessMap.insert( ProcessPair( G4HyperAlpha::Definition(),
                                   theHyperAlphaInelasticProcess ) );
  G4HadronicProcess* theAntiHyperAlphaInelasticProcess =
    new G4HadronInelasticProcess( "anti_hyperalphaInelastic", G4AntiHyperAlpha::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiHyperAlpha::Definition(),
                                   theAntiHyperAlphaInelasticProcess ) );
  G4HadronicProcess* theHyperH4InelasticProcess =
    new G4HadronInelasticProcess( "hyperH4Inelastic", G4HyperH4::Definition() );
  fProcessMap.insert( ProcessPair( G4HyperH4::Definition(),
                                   theHyperH4InelasticProcess ) );
  G4HadronicProcess* theAntiHyperH4InelasticProcess =
    new G4HadronInelasticProcess( "anti_hyperH4Inelastic", G4AntiHyperH4::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiHyperH4::Definition(),
                                   theAntiHyperH4InelasticProcess ) );
  G4HadronicProcess* theDoubleHyperH4InelasticProcess =
    new G4HadronInelasticProcess( "doublehyperH4Inelastic", G4DoubleHyperH4::Definition() );
  fProcessMap.insert( ProcessPair( G4DoubleHyperH4::Definition(),
                                   theDoubleHyperH4InelasticProcess ) );
  G4HadronicProcess* theAntiDoubleHyperH4InelasticProcess =
    new G4HadronInelasticProcess( "anti_doublehyperH4Inelastic",
                                  G4AntiDoubleHyperH4::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiDoubleHyperH4::Definition(),
                                   theAntiDoubleHyperH4InelasticProcess ) );
  G4HadronicProcess* theDoubleHyperDoubleNeutronInelasticProcess =
    new G4HadronInelasticProcess( "doublehyperdoubleneutronInelastic",
                                  G4DoubleHyperDoubleNeutron::Definition() );
  fProcessMap.insert( ProcessPair( G4DoubleHyperDoubleNeutron::Definition(),
                                   theDoubleHyperDoubleNeutronInelasticProcess ) );
  G4HadronicProcess* theAntiDoubleHyperDoubleNeutronInelasticProcess =
    new G4HadronInelasticProcess( "anti_doublehyperdoubleneutronInelastic",
                                  G4AntiDoubleHyperDoubleNeutron::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiDoubleHyperDoubleNeutron::Definition(),
                                   theAntiDoubleHyperDoubleNeutronInelasticProcess ) );
  G4HadronicProcess* theHyperHe5InelasticProcess =
    new G4HadronInelasticProcess( "hyperHe5Inelastic", G4HyperHe5::Definition() );
  fProcessMap.insert( ProcessPair( G4HyperHe5::Definition(),
                                   theHyperHe5InelasticProcess ) );
  G4HadronicProcess* theAntiHyperHe5InelasticProcess =
    new G4HadronInelasticProcess( "anti_hyperHe5Inelastic", G4AntiHyperHe5::Definition() );
  fProcessMap.insert( ProcessPair( G4AntiHyperHe5::Definition(),
                                   theAntiHyperHe5InelasticProcess ) );
  
  // Add the cross sections to the corresponding hadronic processes
  thePionMinusInelasticProcess->AddDataSet( thePionMinusXSdata );
  thePionPlusInelasticProcess->AddDataSet( thePionPlusXSdata );
  theKaonMinusInelasticProcess->AddDataSet( theKaonXSdata );
  theKaonPlusInelasticProcess->AddDataSet( theKaonXSdata );
  theKaonZeroLInelasticProcess->AddDataSet( theKaonXSdata );
  theKaonZeroSInelasticProcess->AddDataSet( theKaonXSdata );
  theProtonInelasticProcess->AddDataSet( theProtonXSdata );
  theNeutronInelasticProcess->AddDataSet( theNeutronXSdata );
  theDeuteronInelasticProcess->AddDataSet( theNuclNuclXSdata );
  theTritonInelasticProcess->AddDataSet( theNuclNuclXSdata );
  theHe3InelasticProcess->AddDataSet( theNuclNuclXSdata );
  theAlphaInelasticProcess->AddDataSet( theNuclNuclXSdata );
  theIonInelasticProcess->AddDataSet( theNuclNuclXSdata );
  theLambdaInelasticProcess->AddDataSet( theHyperonsXSdata );
  theSigmaMinusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theSigmaPlusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theXiMinusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theXiZeroInelasticProcess->AddDataSet( theHyperonsXSdata );
  theOmegaMinusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theAntiProtonInelasticProcess->AddDataSet( theAntibaryonsXSdata );
  theAntiNeutronInelasticProcess->AddDataSet( theAntibaryonsXSdata );
  theAntiDeuteronInelasticProcess->AddDataSet( theAntibaryonsXSdata );
  theAntiTritonInelasticProcess->AddDataSet( theAntibaryonsXSdata );
  theAntiHe3InelasticProcess->AddDataSet( theAntibaryonsXSdata );
  theAntiAlphaInelasticProcess->AddDataSet( theHyperonsXSdata );
  theAntiLambdaInelasticProcess->AddDataSet( theHyperonsXSdata );
  theAntiSigmaMinusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theAntiSigmaPlusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theAntiXiMinusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theAntiXiZeroInelasticProcess->AddDataSet( theHyperonsXSdata );
  theAntiOmegaMinusInelasticProcess->AddDataSet( theHyperonsXSdata );

  theDPlusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theDMinusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theDZeroInelasticProcess->AddDataSet( theHyperonsXSdata );
  theAntiDZeroInelasticProcess->AddDataSet( theHyperonsXSdata );
  theDsPlusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theDsMinusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theBPlusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theBMinusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theBZeroInelasticProcess->AddDataSet( theHyperonsXSdata );
  theAntiBZeroInelasticProcess->AddDataSet( theHyperonsXSdata );
  theBsZeroInelasticProcess->AddDataSet( theHyperonsXSdata );
  theAntiBsZeroInelasticProcess->AddDataSet( theHyperonsXSdata );
  theBcPlusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theBcMinusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theLambdacPlusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theAntiLambdacPlusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theXicPlusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theAntiXicPlusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theXicZeroInelasticProcess->AddDataSet( theHyperonsXSdata );
  theAntiXicZeroInelasticProcess->AddDataSet( theHyperonsXSdata );
  theOmegacZeroInelasticProcess->AddDataSet( theHyperonsXSdata );
  theAntiOmegacZeroInelasticProcess->AddDataSet( theHyperonsXSdata );
  theLambdabInelasticProcess->AddDataSet( theHyperonsXSdata );
  theAntiLambdabInelasticProcess->AddDataSet( theHyperonsXSdata );
  theXibZeroInelasticProcess->AddDataSet( theHyperonsXSdata );
  theAntiXibZeroInelasticProcess->AddDataSet( theHyperonsXSdata );
  theXibMinusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theAntiXibMinusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theOmegabMinusInelasticProcess->AddDataSet( theHyperonsXSdata );
  theAntiOmegabMinusInelasticProcess->AddDataSet( theHyperonsXSdata );

  theHyperTritonInelasticProcess->AddDataSet( theNuclNuclXSdata );
  theAntiHyperTritonInelasticProcess->AddDataSet( theAntibaryonsXSdata );
  theHyperAlphaInelasticProcess->AddDataSet( theNuclNuclXSdata );
  theAntiHyperAlphaInelasticProcess->AddDataSet( theAntibaryonsXSdata );
  theHyperH4InelasticProcess->AddDataSet( theNuclNuclXSdata );
  theAntiHyperH4InelasticProcess->AddDataSet( theAntibaryonsXSdata );
  theDoubleHyperH4InelasticProcess->AddDataSet( theNuclNuclXSdata );
  theAntiDoubleHyperH4InelasticProcess->AddDataSet( theAntibaryonsXSdata );
  theDoubleHyperDoubleNeutronInelasticProcess->AddDataSet( theNuclNuclXSdata );
  theAntiDoubleHyperDoubleNeutronInelasticProcess->AddDataSet( theAntibaryonsXSdata );
  theHyperHe5InelasticProcess->AddDataSet( theNuclNuclXSdata );
  theAntiHyperHe5InelasticProcess->AddDataSet( theAntibaryonsXSdata );

  // Register the proper hadronic model(s) to the corresponding hadronic processes.
  // Note: hadronic models ("BERT", "BIC", "IonBIC", "INCL", "FTFP", "QGSP") are
  //       used for the hadrons and energies they are applicable
  //       (exception for INCL, which in recent versions of Geant4 can handle
  //        more hadron types and higher energies than considered here).
  //       For "physics-list proxies" ("FTFP_BERT", "FTFP_BERT_ATL", "QGSP_BERT",
  //       "QGSP_BIC", "FTFP_INCLXX"), all hadron types and all energies are covered
  //       by combining different hadronic models - similarly (but not identically)
  //       to the corresponding physics lists.
  if ( fPhysicsCase == "BIC"  ||
       fPhysicsCase == "QGSP_BIC" ) {
    // The BIC model is applicable to nucleons and pions,
    // whereas in the physics list QGSP_BIC it is used only for nucleons
    fPhysicsCaseIsSupported = true;
    theProtonInelasticProcess->RegisterMe( theBICmodel );
    theNeutronInelasticProcess->RegisterMe( theBICmodel );
    if ( fPhysicsCase == "BIC" ) {
      thePionMinusInelasticProcess->RegisterMe( theBICmodel );
      thePionPlusInelasticProcess->RegisterMe( theBICmodel );
    } else {
      thePionMinusInelasticProcess->RegisterMe( theBERTmodel );
      thePionPlusInelasticProcess->RegisterMe( theBERTmodel );
    }
  } else if ( fPhysicsCase == "INCL"  ||
              fPhysicsCase == "FTFP_INCLXX" ) {
    // We consider here for simplicity only nucleons and pions
    // (although recent versions of INCL can handle others particles as well)
    fPhysicsCaseIsSupported = true;
    thePionMinusInelasticProcess->RegisterMe( theINCLmodel );
    thePionPlusInelasticProcess->RegisterMe( theINCLmodel );
    theProtonInelasticProcess->RegisterMe( theINCLmodel );
    theNeutronInelasticProcess->RegisterMe( theINCLmodel );
  }
  if ( fPhysicsCase == "IonBIC"         ||
       fPhysicsCase == "FTFP_BERT_ATL"  ||
       fPhysicsCase == "FTFP_BERT"      ||
       fPhysicsCase == "FTFP_INCLXX"    ||
       fPhysicsCase == "QGSP_BERT"      ||
       fPhysicsCase == "QGSP_BIC"  ) {
    // The Binary Light Ion model is used for ions in all physics lists
    fPhysicsCaseIsSupported = true;
    theDeuteronInelasticProcess->RegisterMe( theIonBICmodel );
    theTritonInelasticProcess->RegisterMe( theIonBICmodel );
    theHe3InelasticProcess->RegisterMe( theIonBICmodel );
    theAlphaInelasticProcess->RegisterMe( theIonBICmodel );
    theIonInelasticProcess->RegisterMe( theIonBICmodel );
  }
  if ( fPhysicsCase == "QGSP"       ||
       fPhysicsCase == "QGSP_BERT"  ||
       fPhysicsCase == "QGSP_BIC" ) {
    fPhysicsCaseIsSupported = true;
    thePionMinusInelasticProcess->RegisterMe( theQGSPmodel );
    thePionPlusInelasticProcess->RegisterMe( theQGSPmodel );
    theKaonMinusInelasticProcess->RegisterMe( theQGSPmodel );
    theKaonPlusInelasticProcess->RegisterMe( theQGSPmodel );
    theKaonZeroLInelasticProcess->RegisterMe( theQGSPmodel );
    theKaonZeroSInelasticProcess->RegisterMe( theQGSPmodel );
    theProtonInelasticProcess->RegisterMe( theQGSPmodel );
    theNeutronInelasticProcess->RegisterMe( theQGSPmodel );
    theLambdaInelasticProcess->RegisterMe( theQGSPmodel );
    theSigmaMinusInelasticProcess->RegisterMe( theQGSPmodel );
    theSigmaPlusInelasticProcess->RegisterMe( theQGSPmodel );
    theXiMinusInelasticProcess->RegisterMe( theQGSPmodel );
    theXiZeroInelasticProcess->RegisterMe( theQGSPmodel );
    theOmegaMinusInelasticProcess->RegisterMe( theQGSPmodel );
    theAntiProtonInelasticProcess->RegisterMe( theQGSPmodel );
    theAntiNeutronInelasticProcess->RegisterMe( theQGSPmodel );
    theAntiLambdaInelasticProcess->RegisterMe( theQGSPmodel );
    theAntiSigmaMinusInelasticProcess->RegisterMe( theQGSPmodel );
    theAntiSigmaPlusInelasticProcess->RegisterMe( theQGSPmodel );
    theAntiXiMinusInelasticProcess->RegisterMe( theQGSPmodel );
    theAntiXiZeroInelasticProcess->RegisterMe( theQGSPmodel );
    theAntiOmegaMinusInelasticProcess->RegisterMe( theQGSPmodel );
    theDPlusInelasticProcess->RegisterMe( theQGSPmodel );
    theDMinusInelasticProcess->RegisterMe( theQGSPmodel );
    theDZeroInelasticProcess->RegisterMe( theQGSPmodel );
    theAntiDZeroInelasticProcess->RegisterMe( theQGSPmodel );
    theDsPlusInelasticProcess->RegisterMe( theQGSPmodel );
    theDsMinusInelasticProcess->RegisterMe( theQGSPmodel );
    theBPlusInelasticProcess->RegisterMe( theQGSPmodel );
    theBMinusInelasticProcess->RegisterMe( theQGSPmodel );
    theBZeroInelasticProcess->RegisterMe( theQGSPmodel );
    theAntiBZeroInelasticProcess->RegisterMe( theQGSPmodel );
    theBsZeroInelasticProcess->RegisterMe( theQGSPmodel );
    theAntiBsZeroInelasticProcess->RegisterMe( theQGSPmodel );
    theBcPlusInelasticProcess->RegisterMe( theQGSPmodel );
    theBcMinusInelasticProcess->RegisterMe( theQGSPmodel );
    theLambdacPlusInelasticProcess->RegisterMe( theQGSPmodel );
    theAntiLambdacPlusInelasticProcess->RegisterMe( theQGSPmodel );
    theXicPlusInelasticProcess->RegisterMe( theQGSPmodel );
    theAntiXicPlusInelasticProcess->RegisterMe( theQGSPmodel );
    theXicZeroInelasticProcess->RegisterMe( theQGSPmodel );
    theAntiXicZeroInelasticProcess->RegisterMe( theQGSPmodel );
    theOmegacZeroInelasticProcess->RegisterMe( theQGSPmodel );
    theAntiOmegacZeroInelasticProcess->RegisterMe( theQGSPmodel );
    theLambdabInelasticProcess->RegisterMe( theQGSPmodel );
    theAntiLambdabInelasticProcess->RegisterMe( theQGSPmodel );
    theXibZeroInelasticProcess->RegisterMe( theQGSPmodel );
    theAntiXibZeroInelasticProcess->RegisterMe( theQGSPmodel );
    theXibMinusInelasticProcess->RegisterMe( theQGSPmodel );
    theAntiXibMinusInelasticProcess->RegisterMe( theQGSPmodel );
    theOmegabMinusInelasticProcess->RegisterMe( theQGSPmodel );
    theAntiOmegabMinusInelasticProcess->RegisterMe( theQGSPmodel );
  }
  if ( fPhysicsCase == "BERT"           ||
       fPhysicsCase == "FTFP_BERT_ATL"  ||
       fPhysicsCase == "FTFP_BERT"      ||
       fPhysicsCase == "QGSP_BERT" ) {
    // The BERT model is used for pions and nucleons in all Bertini-based physics lists
    fPhysicsCaseIsSupported = true;
    thePionMinusInelasticProcess->RegisterMe( theBERTmodel );
    thePionPlusInelasticProcess->RegisterMe( theBERTmodel );
    theProtonInelasticProcess->RegisterMe( theBERTmodel );
    theNeutronInelasticProcess->RegisterMe( theBERTmodel );
  }
  if ( fPhysicsCase == "BERT"           ||
       fPhysicsCase == "FTFP_BERT_ATL"  ||
       fPhysicsCase == "FTFP_BERT"      ||
       fPhysicsCase == "FTFP_INCLXX"    ||
       fPhysicsCase == "QGSP_BERT"      ||
       fPhysicsCase == "QGSP_BIC" ) {
    // The BERT model is used for kaons and hyperons in all physics lists
    fPhysicsCaseIsSupported = true;
    theKaonMinusInelasticProcess->RegisterMe( theBERTmodel );
    theKaonPlusInelasticProcess->RegisterMe( theBERTmodel );
    theKaonZeroLInelasticProcess->RegisterMe( theBERTmodel );
    theKaonZeroSInelasticProcess->RegisterMe( theBERTmodel );
    theLambdaInelasticProcess->RegisterMe( theBERTmodel );
    theSigmaMinusInelasticProcess->RegisterMe( theBERTmodel );
    theSigmaPlusInelasticProcess->RegisterMe( theBERTmodel );
    theXiMinusInelasticProcess->RegisterMe( theBERTmodel );
    theXiZeroInelasticProcess->RegisterMe( theBERTmodel );
    theOmegaMinusInelasticProcess->RegisterMe( theBERTmodel );
  }
  if ( fPhysicsCase == "FTFP"           ||
       fPhysicsCase == "FTFP_BERT_ATL"  ||
       fPhysicsCase == "FTFP_BERT"      ||
       fPhysicsCase == "FTFP_INCLXX"    ||
       fPhysicsCase == "QGSP_BERT"      ||
       fPhysicsCase == "QGSP_BIC" ) {
    // The FTFP model is applied for all hadrons, but in different energy intervals according
    // whether it is consider as a stand-alone hadronic model, or within physics lists
    fPhysicsCaseIsSupported = true;
    theAntiDeuteronInelasticProcess->RegisterMe( theFTFPmodel );
    theAntiTritonInelasticProcess->RegisterMe( theFTFPmodel );
    theAntiHe3InelasticProcess->RegisterMe( theFTFPmodel );
    theAntiAlphaInelasticProcess->RegisterMe( theFTFPmodel );
    G4TheoFSGenerator* theFTFPmodelToBeUsed = theFTFPmodel_aboveThreshold;
    if ( fPhysicsCase == "FTFP" ) {
      theFTFPmodelToBeUsed = theFTFPmodel;
    } else if ( fPhysicsCase == "QGSP_BERT"  ||  fPhysicsCase == "QGSP_BIC" ) {
      theFTFPmodelToBeUsed = theFTFPmodel_constrained;
    }
    thePionMinusInelasticProcess->RegisterMe( theFTFPmodelToBeUsed );
    thePionPlusInelasticProcess->RegisterMe( theFTFPmodelToBeUsed );
    theKaonMinusInelasticProcess->RegisterMe( theFTFPmodelToBeUsed );
    theKaonPlusInelasticProcess->RegisterMe( theFTFPmodelToBeUsed );
    theKaonZeroLInelasticProcess->RegisterMe( theFTFPmodelToBeUsed );
    theKaonZeroSInelasticProcess->RegisterMe( theFTFPmodelToBeUsed );
    theProtonInelasticProcess->RegisterMe( theFTFPmodelToBeUsed );
    theAntiProtonInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theNeutronInelasticProcess->RegisterMe( theFTFPmodelToBeUsed );
    theAntiNeutronInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theLambdaInelasticProcess->RegisterMe( theFTFPmodelToBeUsed );
    theAntiLambdaInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theSigmaMinusInelasticProcess->RegisterMe( theFTFPmodelToBeUsed );
    theAntiSigmaMinusInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theSigmaPlusInelasticProcess->RegisterMe( theFTFPmodelToBeUsed );
    theAntiSigmaPlusInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theXiMinusInelasticProcess->RegisterMe( theFTFPmodelToBeUsed );
    theAntiXiMinusInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theXiZeroInelasticProcess->RegisterMe( theFTFPmodelToBeUsed );
    theAntiXiZeroInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theOmegaMinusInelasticProcess->RegisterMe( theFTFPmodelToBeUsed );
    theAntiOmegaMinusInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theDPlusInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theDMinusInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theDZeroInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theAntiDZeroInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theDsPlusInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theDsMinusInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theBPlusInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theBMinusInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theBZeroInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theAntiBZeroInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theBsZeroInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theAntiBsZeroInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theBcPlusInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theBcMinusInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theLambdacPlusInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theAntiLambdacPlusInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theXicPlusInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theAntiXicPlusInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theXicZeroInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theAntiXicZeroInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theOmegacZeroInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theAntiOmegacZeroInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theLambdabInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theAntiLambdabInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theXibZeroInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theAntiXibZeroInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theXibMinusInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theAntiXibMinusInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theOmegabMinusInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theAntiOmegabMinusInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    theFTFPmodelToBeUsed = theFTFPmodel_aboveThreshold;
    if ( fPhysicsCase == "FTFP" ) theFTFPmodelToBeUsed = theFTFPmodel;
    theDeuteronInelasticProcess->RegisterMe( theFTFPmodelToBeUsed );
    theTritonInelasticProcess->RegisterMe( theFTFPmodelToBeUsed );
    theHe3InelasticProcess->RegisterMe( theFTFPmodelToBeUsed );
    theAlphaInelasticProcess->RegisterMe( theFTFPmodelToBeUsed );
    theIonInelasticProcess->RegisterMe( theFTFPmodelToBeUsed );
  }

  if ( G4HadronicParameters::Instance()->EnableHyperNuclei() ) {
    // Only FTFP and INCL can handle the nuclear interactions
    // of light hypernuclei, and only FTFP is capable of handling
    // the nuclear interactions of light anti-hypernuclei.
    if ( fPhysicsCase == "FTFP_BERT"    ||
         fPhysicsCase == "FTFP_INCLXX"  ||
         fPhysicsCase == "FTFP" ) {
      fPhysicsCaseIsSupported = true;
      if ( fPhysicsCase == "FTFP_INCLXX" ) {
        theHyperTritonInelasticProcess->RegisterMe( theFTFPmodel );
        theHyperAlphaInelasticProcess->RegisterMe( theFTFPmodel );
        theHyperH4InelasticProcess->RegisterMe( theFTFPmodel );
        theDoubleHyperH4InelasticProcess->RegisterMe( theFTFPmodel );
        theDoubleHyperDoubleNeutronInelasticProcess->RegisterMe( theFTFPmodel );
        theHyperHe5InelasticProcess->RegisterMe( theFTFPmodel );
      } else {
        theHyperTritonInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
        theHyperAlphaInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
        theHyperH4InelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
        theDoubleHyperH4InelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
        theDoubleHyperDoubleNeutronInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
        theHyperHe5InelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
      }
      theAntiHyperTritonInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
      theAntiHyperAlphaInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
      theAntiHyperH4InelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
      theAntiDoubleHyperH4InelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
      theAntiDoubleHyperDoubleNeutronInelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
      theAntiHyperHe5InelasticProcess->RegisterMe( theFTFPmodel_belowThreshold );
    }
    if ( fPhysicsCase == "FTFP_INCLXX"  ||
         fPhysicsCase == "INCL" ) {
      fPhysicsCaseIsSupported = true;
      theHyperTritonInelasticProcess->RegisterMe( theINCLmodel );
      theHyperAlphaInelasticProcess->RegisterMe( theINCLmodel );
      theHyperH4InelasticProcess->RegisterMe( theINCLmodel );
      theDoubleHyperH4InelasticProcess->RegisterMe( theINCLmodel );
      theDoubleHyperDoubleNeutronInelasticProcess->RegisterMe( theINCLmodel );
      theHyperHe5InelasticProcess->RegisterMe( theINCLmodel );
    }
  }

  if ( ! fPhysicsCaseIsSupported ) {
    G4cerr << "ERROR: Not supported final-state hadronic inelastic physics case !"
           << fPhysicsCase << G4endl
           << "\t Re-try by choosing one of the following:" << G4endl
           << "\t - Hadronic models : BERT, BIC, IonBIC, INCL, FTFP, QGSP" << G4endl
           << "\t - \"Physics-list proxies\" : FTFP_BERT (default), FTFP_BERT_ATL, \
                                               QGSP_BERT, QGSP_BIC, FTFP_INCLXX"
           << G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HadronicGenerator::~HadronicGenerator() {
  fPartTable->DeleteAllParticles();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool HadronicGenerator::IsApplicable( const G4String &nameProjectile,
                                        const G4double projectileEnergy ) const {
  G4ParticleDefinition* projectileDefinition = fPartTable->FindParticle( nameProjectile );
  return IsApplicable( projectileDefinition, projectileEnergy );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool HadronicGenerator::IsApplicable( G4ParticleDefinition* projectileDefinition,
                                        const G4double projectileEnergy ) const {
  if ( projectileDefinition == nullptr ) return false;
  G4bool isApplicable = true;
  // No restrictions for "physics list proxies" because they cover all hadron types and energies.
  // For the individual models, instead, we need to consider their limitations. 
  if ( fPhysicsCase == "BERT" ) {
    // We consider BERT model below 15 GeV
    if ( ( ( projectileDefinition != G4PionMinus::Definition()  ) &&
           ( projectileDefinition != G4PionPlus::Definition()   ) &&
           ( projectileDefinition != G4Proton::Definition()     ) &&
           ( projectileDefinition != G4Neutron::Definition()    ) &&
           ( projectileDefinition != G4Lambda::Definition()     ) &&
           ( projectileDefinition != G4SigmaMinus::Definition() ) &&
           ( projectileDefinition != G4SigmaPlus::Definition()  ) &&
           ( projectileDefinition != G4XiMinus::Definition()    ) &&
           ( projectileDefinition != G4XiZero::Definition()     ) &&
           ( projectileDefinition != G4OmegaMinus::Definition() ) ) ||
         ( projectileEnergy > 15.0*CLHEP::GeV ) ) {
      isApplicable = false;
    }
  } else if ( fPhysicsCase == "QGSP" ) {
    // We consider QGSP above 2 GeV and not for ions or anti-ions
    if ( projectileEnergy < 2.0*CLHEP::GeV                    ||
         projectileDefinition == G4Deuteron::Definition()     ||
         projectileDefinition == G4Triton::Definition()       ||
         projectileDefinition == G4He3::Definition()          ||
         projectileDefinition == G4Alpha::Definition()        ||
         projectileDefinition == G4GenericIon::Definition()   ||
         projectileDefinition == G4AntiDeuteron::Definition() ||
         projectileDefinition == G4AntiTriton::Definition()   ||
         projectileDefinition == G4AntiHe3::Definition()      ||
         projectileDefinition == G4AntiAlpha::Definition() ) {
      isApplicable = false;
    }
  } else if ( fPhysicsCase == "BIC"  ||  fPhysicsCase == "INCL" ) {
    // We consider BIC and INCL models only for pions and nucleons below 10 GeV
    // (although in recent versions INCL is capable of handling more hadrons
    // and up to higher energies)
    if ( ( ( projectileDefinition != G4PionMinus::Definition() ) &&
           ( projectileDefinition != G4PionPlus::Definition()  ) &&
           ( projectileDefinition != G4Proton::Definition()    ) &&
           ( projectileDefinition != G4Neutron::Definition()   ) ) ||
         ( projectileEnergy > 10.0*CLHEP::GeV ) ) {
      isApplicable = false;
    }
  } else if ( fPhysicsCase == "IonBIC" ) {
    // We consider IonBIC models only for deuteron, triton, He3, alpha
    // with energies below 10 GeV / nucleon
    if ( ! ( ( projectileDefinition == G4Deuteron::Definition()  &&
               projectileEnergy < 2*10.0*CLHEP::GeV                 ) ||
             ( projectileDefinition == G4Triton::Definition()    &&
               projectileEnergy < 3*10.0*CLHEP::GeV                 ) ||
             ( projectileDefinition == G4He3::Definition()       &&
               projectileEnergy < 3*10.0*CLHEP::GeV                 ) ||
             ( projectileDefinition == G4Alpha::Definition()     &&
               projectileEnergy < 4*10.0*CLHEP::GeV                 ) ) ) {
      isApplicable = false;
    }
  } 
  return isApplicable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* HadronicGenerator::
GenerateInteraction( const G4String &nameProjectile, const G4double projectileEnergy,
                     const G4ThreeVector &projectileDirection, G4Material* targetMaterial ) {
  G4ParticleDefinition* projectileDefinition = fPartTable->FindParticle( nameProjectile );
  return GenerateInteraction( projectileDefinition, projectileEnergy,
                              projectileDirection, targetMaterial );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* HadronicGenerator::
GenerateInteraction( G4ParticleDefinition* projectileDefinition, const G4double projectileEnergy,
                     const G4ThreeVector &projectileDirection, G4Material* targetMaterial ) {
  // This is the most important method of the HadronicGenerator class:
  // the method performs the specified hadronic interaction
  // (by invoking the "PostStepDoIt" method of the corresponding hadronic process)
  // and returns the final state, i.e. the secondaries produced by the collision.
  // It is a relatively short method because the heavy load of setting up all
  // possible hadronic processes - with their hadronic models, transition regions,
  // and cross sections (the latter is needed for sampling the target nucleus from
  // the target material) - was already done by the constructor of the class.
  G4VParticleChange* aChange = nullptr;

  if ( projectileDefinition == nullptr ) {
    G4cerr << "ERROR: projectileDefinition is NULL !" << G4endl;
    return aChange;
  }

  // Debugging print-out
  //G4cout << "\t" << projectileDefinition->GetParticleName()
  //       << "\t" << projectileEnergy/CLHEP::GeV
  //       << " GeV \t" << projectileDirection
  //       << "\t" << ( targetMaterial ? targetMaterial->GetName() : "NULL" );
  if ( ! IsApplicable( projectileDefinition, projectileEnergy ) ) {
    //G4cout << " -> NOT applicable !" ; //<< G4endl;  // Debugging print-out
    return aChange;
  }
  //G4cout << G4endl;

  // Check Geant4 state (not strictly needed)
  //if ( ! G4StateManager::GetStateManager()->SetNewState( G4State_PreInit ) ) {
  //  G4cerr << "ERROR: No possible to set G4State_PreInit !" << G4endl;
  //  return aChange;
  //}

  // Geometry definition (not strictly needed)
  //const G4double dimX = 1.0*mm;
  //const G4double dimY = 1.0*mm;
  //const G4double dimZ = 1.0*mm;
  //G4Box* sFrame = new G4Box( "Box", dimX, dimY, dimZ );
  //G4LogicalVolume* lFrame = new G4LogicalVolume( sFrame, targetMaterial, "Box", 0, 0, 0 );
  //G4PVPlacement* pFrame = new G4PVPlacement( 0, G4ThreeVector(), "Box", lFrame, 0, false, 0 );
  //G4TransportationManager::GetTransportationManager()->SetWorldForTracking( pFrame );

  // Projectile track & step
  G4DynamicParticle dParticle( projectileDefinition, projectileDirection, projectileEnergy );
  const G4double aTime = 0.0;
  const G4ThreeVector aPosition = G4ThreeVector( 0.0, 0.0, 0.0 );
  G4Track* gTrack = new G4Track( &dParticle, aTime, aPosition );
  G4TouchableHandle fpTouchable( new G4TouchableHistory );  // Not strictly needed
  gTrack->SetTouchableHandle( fpTouchable );                // Not strictly needed
  G4Step* step = new G4Step;
  step->SetTrack( gTrack );
  gTrack->SetStep( step );
  G4StepPoint* aPoint = new G4StepPoint;
  aPoint->SetPosition( aPosition );
  aPoint->SetMaterial( targetMaterial );
  step->SetPreStepPoint( aPoint );
  dParticle.SetKineticEnergy( projectileEnergy );
  gTrack->SetStep( step );
  gTrack->SetKineticEnergy( projectileEnergy );

  // Change Geant4 state: from "PreInit" to "Idle" (not strictly needed)
  //if ( ! G4StateManager::GetStateManager()->SetNewState( G4State_Idle ) ) {
  //  G4cerr << "ERROR: No possible to set G4State_Idle !" << G4endl;
  //  return aChange;
  //}

  // Finally, the hadronic interaction: hadron projectile and ion projectile
  // need to be treated slightly differently
  G4HadronicProcess* theProcess = nullptr;
  G4ParticleDefinition* theProjectileDef = nullptr;
  if ( projectileDefinition->IsGeneralIon() ) {
    theProjectileDef = G4GenericIon::Definition();
  } else {
    theProjectileDef = projectileDefinition;
  }
  auto mapIndex = fProcessMap.find( theProjectileDef );
  if ( mapIndex != fProcessMap.end() ) theProcess = mapIndex->second;
  if ( theProcess != nullptr ) {
    aChange = theProcess->PostStepDoIt( *gTrack, *step );
    //**************************************************
  } else {
    G4cerr << "ERROR: theProcess is nullptr !" << G4endl;
  }
  fLastHadronicProcess = theProcess;
  //delete pFrame;
  //delete lFrame;
  //delete sFrame;

  return aChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double HadronicGenerator::GetImpactParameter() const {
  G4double impactParameter = -999.0 * fermi;
  G4HadronicProcess* hadProcess = GetHadronicProcess();
  G4HadronicInteraction* hadInteraction = GetHadronicInteraction();
  G4HadronicInteraction* wantedHadInteraction =
    const_cast< G4HadronicProcess* >( hadProcess )->GetHadronicModel( "FTFP" );
  if ( hadInteraction != nullptr && hadInteraction == wantedHadInteraction ) {
    // FTFP has handled the inelastic hadronic interaction.
    G4TheoFSGenerator* theoFSGenerator = dynamic_cast< G4TheoFSGenerator* >( hadInteraction );
    if ( theoFSGenerator != nullptr ) {
      const G4FTFModel* ftfModel =
        dynamic_cast< const G4FTFModel* >( theoFSGenerator->GetHighEnergyGenerator() );
      if ( ftfModel != nullptr ) {
        // ftfModel points to the G4FTFModel object instance that handled the
        // inelastic hadronic interaction.
        impactParameter = ftfModel->GetImpactParameter();
        //G4cout << "\t impactParameter = " << impactParameter/fermi << " fm" << G4endl;
      }
    }
  }
  return impactParameter;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int HadronicGenerator::GetNumberOfProjectileSpectatorNucleons() const {
  G4double numProjectileSpectatorNucleons = -999;
  G4HadronicProcess* hadProcess = GetHadronicProcess();
  G4HadronicInteraction* hadInteraction = GetHadronicInteraction();
  G4HadronicInteraction* wantedHadInteraction =
    const_cast< G4HadronicProcess* >( hadProcess )->GetHadronicModel( "FTFP" );
  if ( hadInteraction != nullptr && hadInteraction == wantedHadInteraction ) {
    G4TheoFSGenerator* theoFSGenerator = dynamic_cast< G4TheoFSGenerator* >( hadInteraction );
    if ( theoFSGenerator != nullptr ) {
      const G4FTFModel* ftfModel =
        dynamic_cast< const G4FTFModel* >( theoFSGenerator->GetHighEnergyGenerator() );
      if ( ftfModel != nullptr ) {
        numProjectileSpectatorNucleons = ftfModel->GetNumberOfProjectileSpectatorNucleons();
        //G4cout << "\t numProjectileSpectatorNucleons = " << numProjectileSpectatorNucleons
        //       << G4endl;
      }
    }
  }
  return numProjectileSpectatorNucleons;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int HadronicGenerator::GetNumberOfTargetSpectatorNucleons() const {
  G4double numTargetSpectatorNucleons = -999;
  G4HadronicProcess* hadProcess = GetHadronicProcess();
  G4HadronicInteraction* hadInteraction = GetHadronicInteraction();
  G4HadronicInteraction* wantedHadInteraction =
    const_cast< G4HadronicProcess* >( hadProcess )->GetHadronicModel( "FTFP" );
  if ( hadInteraction != nullptr && hadInteraction == wantedHadInteraction ) {
    G4TheoFSGenerator* theoFSGenerator = dynamic_cast< G4TheoFSGenerator* >( hadInteraction );
    if ( theoFSGenerator != nullptr ) {
      const G4FTFModel* ftfModel =
        dynamic_cast< const G4FTFModel* >( theoFSGenerator->GetHighEnergyGenerator() );
      if ( ftfModel != nullptr ) {
        numTargetSpectatorNucleons = ftfModel->GetNumberOfTargetSpectatorNucleons();
        //G4cout << "\t numTargetSpectatorNucleons = " << numTargetSpectatorNucleons << G4endl;
      }
    }
  }
  return numTargetSpectatorNucleons;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int HadronicGenerator::GetNumberOfNNcollisions() const {
  G4double numNNcollisions = -999;
  G4HadronicProcess* hadProcess = GetHadronicProcess();
  G4HadronicInteraction* hadInteraction = GetHadronicInteraction();
  G4HadronicInteraction* wantedHadInteraction =
    const_cast< G4HadronicProcess* >( hadProcess )->GetHadronicModel( "FTFP" );
  if ( hadInteraction != nullptr && hadInteraction == wantedHadInteraction ) {
    G4TheoFSGenerator* theoFSGenerator = dynamic_cast< G4TheoFSGenerator* >( hadInteraction );
    if ( theoFSGenerator != nullptr ) {
      const G4FTFModel* ftfModel =
        dynamic_cast< const G4FTFModel* >( theoFSGenerator->GetHighEnergyGenerator() );
      if ( ftfModel != nullptr ) {
        numNNcollisions = ftfModel->GetNumberOfNNcollisions();
        //G4cout << "\t numNNcollisions = " << numNNcollisions << G4endl;
      }
    }
  }
  return numNNcollisions;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
