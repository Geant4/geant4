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
// G4PhysicsModelCatalog class implementation
//
// Author: M.Asai (SLAC), 26 September 2013
//
// Revised in August 2021 by A.Ribon (CERN).
// --------------------------------------------------------------------------

#include "G4PhysicsModelCatalog.hh"
#include "G4Threading.hh"

G4bool G4PhysicsModelCatalog::isInitialized = false;
std::vector< G4int >*    G4PhysicsModelCatalog::theVectorOfModelIDs   = nullptr;
std::vector< G4String >* G4PhysicsModelCatalog::theVectorOfModelNames = nullptr;

// --------------------------------------------------------------------------
void G4PhysicsModelCatalog::Initialize() {
  if(isInitialized)
  {
    return;
  }
  if ( theVectorOfModelIDs == nullptr  &&  theVectorOfModelNames == nullptr ) {
    static std::vector< G4int > aVectorOfInts;
    theVectorOfModelIDs = &aVectorOfInts;
    static std::vector< G4String > aVectorOfStrings;
    theVectorOfModelNames = &aVectorOfStrings;

    // NOTE:
    // The goal is that, starting from Geant4 version 11.0, all the three
    // identifiers (modelID, index, name) remain the same regardless of the
    // physics list, application, and version of Geant4.
    // Therefore, after Geant4 11.0, you can only add an entry - e.g. when
    // a new model is added, or when a pre-existing model not yet present
    // in this catalogue is included - at the bottom of this method
    // (rather than inserting it in the middle), and do NOT delete anything.
    //
    // For the modelID values, these are the considerations and choices made:
    // - Values below 1000 are excluded because reserved to process and
    //   sub-process ID values.
    // - Whenever resonable, modelID values should have free spaces between
    //   them, to allow for eventual, future variants - e.g. different
    //   tunings of the same model - to keep modelID values close to the
    //   original model.
    // - modelID values are between 10'000 and 39'999, with the following
    //   subdivision in 3 categories (identified by the most significant
    //   digit):
    //   *  Values between 10'000 and 19'999 are for EM models
    //   *  Values between 20'000 and 29'999 are for HAD models
    //   *  Values between 30'000 and 39'999 are for all the other models
    //      (i.e. included neither in EM models nor in HAD models).
    //   Note that larger values of modelID are neither more difficult to
    //   handle nor less computing performant with respect to smaller values
    //   (we remind that, for plotting, the index of the model, rather than
    //    its modelID, should be conveniently used, whereas for all the rest
    //    the modelID is recommended).
    
    // =======================================================================
    // ================= 1st EM MODELS : from 10'000 to 19'999 ===============
    // =======================================================================

    InsertModel( 10000, "model_EM" );
    
    // e- production
    InsertModel( 10010, "model_DeltaElectron" );
    InsertModel( 10011, "model_DeltaEBelowCut" );
    InsertModel( 10012, "model_PhotoElectron" );
    InsertModel( 10013, "model_ComptonElectron" );
    InsertModel( 10014, "model_TripletElectron" );

    // gamma production
    InsertModel( 10020, "model_Bremsstrahlung" );
    InsertModel( 10021, "model_SplitBremsstrahlung" );
    InsertModel( 10022, "model_ComptonGamma" );
    InsertModel( 10023, "model_Annihilation" );
    InsertModel( 10024, "model_TripletGamma" );
    InsertModel( 10025, "model_GammaGammaEntanglement" );

    // e+e- pair production
    InsertModel( 10030, "model_EplusEminisPair" );

    // atomic de-excitation
    InsertModel( 10040, "model_Fluorescence" );
    InsertModel( 10041, "model_gammaPIXE" );
    InsertModel( 10050, "model_AugerElectron" );
    InsertModel( 10051, "model_ePIXE" );

    // recoil ions
    InsertModel( 10060, "model_IonRecoil" );

    // DNA models
    InsertModel( 11000, "model_DNA" );
    InsertModel( 11001, "model_Ritchie1994eSolvation" );
    InsertModel( 11002, "model_Terrisol1990eSolvation" );
    InsertModel( 11003, "model_Meesungnoen2002eSolvation" );
    InsertModel( 11004, "model_Kreipl2009eSolvation" );
    InsertModel( 11005, "model_MeesungnoenSolid2002eSolvation" );

    // =======================================================================
    // ============= 2nd HADRONIC MODELS : from 20'000 to 29'999 =============
    // =======================================================================

    // ----------------------------------------------------------
    // --- Gamma- , Lepto- , Neutrino-nuclear : 20'000-20'999 ---
    // ----------------------------------------------------------
    //    -  20'000 - 20'099 : Electromagnetic dissociation models
    //    -  20'100 - 20'199 : Gamma-nuclear models
    //    -  20'200 - 20-299 : Electron/positron-nuclear models
    //    -  20'300 - 20'399 : Muon-nuclear models
    //    -  20'400 - 20'499 : Tau-nuclear models
    //    -  20'500 - 20'999 : Neutrino models...
    //    ...

    // --- EM dissociation models: 20'000 - 20'099 ---

    // Class G4EMDissociation
    InsertModel( 20000, "model_projectileEMDissociation" );
    InsertModel( 20001, "model_targetEMDissociation" );

    // --- Gamma-nuclear models: 20'100 - 20'199 ---

    // Class G4LENDorBERTModel
    InsertModel( 20100, "model_LENDorBERTModel" );

    // Class G4LowEGammaNuclearModel
    InsertModel( 20150, "model_GammaNPreco" );

    // --- Charged-lepton - nuclear models: 20'200 - 20'499 ---
    
    // Class G4ElectroVDNuclearModel
    InsertModel( 20200, "model_G4ElectroVDNuclearModel" );
    
    // Class G4MuonVDNuclearModel
    InsertModel( 20300, "model_G4MuonVDNuclearModel" );

    // --- Neutrino models: 20'500 - 20'999 ---
    
    // Class G4NeutrinoElectronCcModel
    InsertModel( 20510, "model_nu-e-inelastic" );

    // Class G4NeutrinoNucleusModel
    InsertModel( 20520, "model_neutrino-nucleus" );
    
    // The following classes derives from G4NeutrinoNucleusModel
    // Class G4ANuElNucleusCcModel
    InsertModel( 20530, "model_ANuElNuclCcModel" );
    // Class G4ANuElNucleusNcModel
    InsertModel( 20540, "model_ANuElNuclNcModel" );
    // Class G4ANuMuNucleusCcModel
    InsertModel( 20550, "model_ANuMuNuclCcModel" );
    // Class G4ANuMuNucleusNcModel
    InsertModel( 20560, "model_ANuMuNuclNcModel" );
    // Class G4NuElNucleusCcModel
    InsertModel( 20570, "model_NuElNuclCcModel" );
    // Class G4NuElNucleusNcModel
    InsertModel( 20580, "model_NuElNuclNcModel" );
    // Class G4NuMuNucleusCcModel
    InsertModel( 20590, "model_NuMuNuclCcModel" );
    // Class G4NuMuNucleusNcModel
    InsertModel( 20600, "model_NuMuNuclNcModel" );
    
    // ------------------------------------------------------------------------------------------
    // --- Elastic, Charge-Exchange, Quasi-Elastic, specialized Diffraction : 21'000 - 21'999 ---
    // ------------------------------------------------------------------------------------------
    //    -  21'000 - 21'199 : Elastic
    //    -  21'200 - 21'299 : Charge-Exchange
    //    -  21'300 - 21'499 : Quasi-Elastic
    //    -  21'500 - 21'999 : specialized Diffraction generators
    // -----------------------------------------------------------

    // --- Elastic models: 21'000 - 21'199 ---
    
    // Class G4HadronElastic
    InsertModel( 21000, "model_hElasticLHEP" );
    // Class G4AntiNuclElastic
    InsertModel( 21010, "model_AntiAElastic" );
    // Class G4ChipsElasticModel
    InsertModel( 21020, "model_hElasticCHIPS" );
    // Class G4DiffuseElastic
    InsertModel( 21030, "model_DiffuseElastic" );
    // Class G4DiffuseElasticV2
    InsertModel( 21040, "model_DiffuseElasticV2" );
    // Class G4NuclNuclDiffuseElastic
    InsertModel( 21050, "model_NNDiffuseElastic" );
    // Class G4ElasticHadrNucleusHE
    InsertModel( 21060, "model_hElasticGlauber" );
    // Class G4hhElastic
    InsertModel( 21070, "model_HadrHadrElastic" );
    // Class G4LowEHadronElastic
    InsertModel( 21080, "model_hLowEElastic" );
    // Class G4LEHadronProtonElastic
    InsertModel( 21090, "model_G4LEHadronProtonElastic" );
    // Class G4LEnp
    InsertModel( 21100, "model_G4LEnp" );
    // Class G4LEpp
    InsertModel( 21110, "model_G4LEpp" );
    // Class G4NeutronElectronElModel
    InsertModel( 21120, "model_n-e-elastic" );
    // Class G4NeutrinoElectronNcModel
    InsertModel( 21130, "model_nu-e-elastic" );
    
    // --- Charge exchange : 21'200 - 21'299 ---
    
    // Class: G4ChargeExchange
    InsertModel( 21200, "model_ChargeExchange" );

    // --- Quasi-Elastic : 21'300 - 21'499 ---
    
    // Class: G4QuasiElasticChannel (which uses Chips's model G4QuasiElRatios)
    InsertModel( 21300, "model_QuasiElastic" );

    // --- Special diffraction generators : 21'500 - 21'999 ---
    
    // Class: G4LMsdGenerator
    InsertModel( 21500, "model_LMsdGenerator" );

    // -----------------------------------------------------------------
    // --- High energy-models (e.g. string models) : 22'000 - 22'999 ---
    // -----------------------------------------------------------------
    //    -  22'000 - 22'099 : G4TheoFSGenerator
    //    -  22'100 - 22'199 : FTF
    //    -  22'200 - 22'299 : QGS
    //    -  ...
    // For gamma         - nuclear : QGS string formation + QGS  string fragmentation + Precompound
    // For e-/e+/mu-/mu+ - nuclear : FTF string formation + Lund string fragmentation + Precompound
    InsertModel( 22000, "model_TheoFSGenerator" );
    // FTF string formation + Lund string fragmentation + Precompound
    InsertModel( 22100, "model_FTFP" );
    // FTF string formation + Lund string fragmentation + Binary Cascade
    InsertModel( 22150, "model_FTFB" );
    // FTF string formation + QGS string fragmentation + Precompound
    InsertModel( 22175, "model_FTFQGSP" );
    // QGS string formation + QGS string fragmentation + Precompound
    InsertModel( 22200, "model_QGSP" );
    // QGS string formation + QGS string fragmentation + Binary Cascade
    InsertModel( 22250, "model_QGSB" );

    // ----------------------------------------------------
    // --- Intermediate energy models : 23'000 - 23'999 ---
    // ----------------------------------------------------
    //    -  23'000 - 23'099 : BERT
    //    -  23'100 - 23'199 : BIC
    //    -  23'200 - 23'299 : INCL
    //    -  23'300 - 23'399 : QMD
    //    ...
    // Class: G4CascadeInterface
    InsertModel( 23000, "model_BertiniCascade" );

    // The names are similar, but not identical, to those before Geant4 11.
    // Class: G4BinaryCascade
    InsertModel( 23100, "model_G4BinaryCascade" );
    // Class: G4BinaryLightIonReaction
    InsertModel( 23110, "model_G4BinaryLightIonReaction" );

    // Class: G4INCLXXInterface
    InsertModel( 23200, "model_INCLXXCascade" );

    // Class: G4QMDReaction
    InsertModel( 23300, "model_QMDModel" );

    // --------------------------------------------------------------
    // --- Pre-equilibrium/De-excitation models : 24'000 - 24'999 ---
    // --------------------------------------------------------------
    //    -  Pre-equilibrium : 24'000 - 24'099
    //       *  24'000 - 24'049 : precompound
    //       *  24'050 - 24'099 : internal BERT pre-equilibrium
    //    -  de-excitation : 24'100 - 24'999
    //       *  24'100 - 24'149 : Evaporation
    //       *  24'150 - 24'199 : Photon evaporation
    //       *  24'200 - 24'299 : GEM evaporation
    //       *  24'300 - 24'349 : Fermi BreakUp
    //       *  24'350 - 24'399 : Multifragmentation
    //       *  24'400 - 24'449 : Ablation
    //       *  24'450 - 24'499 : Fission
    //       *  24'500 - 24'599 : ABLA
    //       *  24'600 - 24'699 : internal BERT de-excitation
    //       ...

    // --- Pre-equilibrium: 24'000 - 24'099 ---

    // Class: G4PreCompoundModel
    InsertModel( 24000, "model_PRECO" );

    // Class: G4LowEIonFragmentation
    InsertModel( 24010, "model_LowEIonPreco" );

    // Class: G4NonEquilibriumEvaporator (i.e. BERT pre-equilibrium)
    InsertModel( 24050, "model_G4NonEquilibriumEvaporator" );

    // --- De-excitation: 24'100 - 24'999 ---

    //     --- Evaporation: 24'100-  24'149 ---
    
    // Class: G4EvaporationChannel
    InsertModel( 24100, "model_G4EvaporationChannel" );

    // Classes that use it: G4NeutronRadCapture and G4ExcitationHandler
    InsertModel( 24110, "model_e-InternalConversion" );
    
    // Class: G4UnstableFragmentBreakUp
    InsertModel( 24120, "model_G4UnstableFragmentBreakUp" );

    //     --- Photon-Evaporation: 24'150 - 24'199 ---
    
    // Class: G4PhotonEvaporation
    InsertModel( 24150, "model_G4PhotonEvaporation" );
    
    // Class: G4NeutronRadCapture
    InsertModel( 24160, "model_nRadCapture" );

    //     --- GEM evaporation : 24'200 - 24'299 ---
    
    // Class: G4GEMChannel
    InsertModel( 24200, "model_G4GEMChannel" );
    
    // Class: G4GEMChannelVI
    InsertModel( 24210, "model_G4GEMChannelVI" );
    
    //     --- Fermi BreakUp : 24'300 - 24'349 ---
    
    // Class: G4FermiBreakUpVI
    InsertModel( 24300, "model_G4FermiBreakUpVI" );
    
    //     --- Multifragmentation : 24'350 - 24'399 ---
    
    // Class: G4StatMF
    InsertModel( 24350, "model_G4StatMF" );

    //     --- Ablation : 24'400 - 24'449 ---
    
    // Class: G4WilsonAblationModel
    InsertModel( 24400, "model_G4WilsonAblationModel" );

    //     --- Fission : 24'450 - 24'499 ---
  
    // Class: G4CompetitiveFission
    InsertModel( 24450, "model_G4CompetitiveFission" );

    //     --- ABLA : 24'500 - 24'599 ---

    // Class: G4AblaInterface
    InsertModel( 24500, "model_ABLA" );
    
    //     --- internal BERT de-excitation : 24'600 - 24'699 ---
    
    // Class: G4EquilibriumEvaporator
    InsertModel( 24600, "model_G4EquilibriumEvaporator" );

    //     --- Other types of de-excitation : 24'700 - 24'999 ---

    //     ...
    
    // ------------------------------------------------
    // --- Low-energy data-driven : 25'000 - 25'999 ---
    // ------------------------------------------------
    //    -  25'000 - 25'199 : ParticleHP
    //    -  25'200 - 25'200 : LEND
    //       ...
    //    -  25'500 - 25'999 : RadioactiveDecay

    // --- ParticleHP : 25'000 - 25'199 ---
    
    // Classes: G4ParticleHPCapture , G4ParticleHPCaptureFS
    InsertModel( 25000, "model_NeutronHPCapture" );

    // Classes: G4ParticleHPElastic , G4ParticleHPElasticFS
    InsertModel( 25010, "model_NeutronHPElastic" );

    // Classes: G4ParticleHPFission , G4ParticleHPFissionFS , G4WendtFissionFragmentGenerator
    InsertModel( 25020, "model_NeutronHPFission" );

    // Inelastic: the following classes inherit from either G4ParticleHPInelasticBaseFS or
    //            G4ParticleHPInelasticCompFS, which in turn inherit from G4ParticleHPFinalState

    // Class G4ParticleHPNInelasticFS
    InsertModel( 25030, "model_G4ParticleHPNInelasticFS_F01" );
    // Class model_G4ParticleHPNXInelasticFS
    InsertModel( 25031, "model_G4ParticleHPNXInelasticFS_F02" );
    // Class G4ParticleHP2NDInelasticFS
    InsertModel( 25032, "model_G4ParticleHP2NDInelasticFS_F03" );
    // Class G4ParticleHP2NInelasticFS
    InsertModel( 25033, "model_G4ParticleHP2NInelasticFS_F04" );
    // Class G4ParticleHP3NInelasticFS
    InsertModel( 25034, "model_G4ParticleHP3NInelasticFS_F05" );
    // Class G4ParticleHPNAInelasticFS
    InsertModel( 25035, "model_G4ParticleHPNAInelasticFS_F06" );
    // Class G4ParticleHPN3AInelasticFS
    InsertModel( 25036, "model_G4ParticleHPN3AInelasticFS_F07" );
    // Class G4ParticleHP2NAInelasticFS
    InsertModel( 25037, "model_G4ParticleHP2NAInelasticFS_F08" );
    // Class G4ParticleHP3NAInelasticFS
    InsertModel( 25038, "model_G4ParticleHP3NAInelasticFS_F09" );
    // Class G4ParticleHPNPInelasticFS
    InsertModel( 25039, "model_G4ParticleHPNPInelasticFS_F10" );
    // Class G4ParticleHPN2AInelasticFS
    InsertModel( 25040, "model_G4ParticleHPN2AInelasticFS_F11" );
    // Clas G4ParticleHP2N2AInelasticFS
    InsertModel( 25041, "model_G4ParticleHP2N2AInelasticFS_F12" );
    // Class G4ParticleHPNDInelasticFS
    InsertModel( 25042, "model_G4ParticleHPNDInelasticFS_F13" );
    // Class G4ParticleHPNTInelasticFS
    InsertModel( 25043, "model_G4ParticleHPNTInelasticFS_F14" );
    // Class G4ParticleHPNHe3InelasticFS
    InsertModel( 25044, "model_G4ParticleHPNHe3InelasticFS_F15" );
    // Class G4ParticleHPND2AInelasticFS
    InsertModel( 25045, "model_G4ParticleHPND2AInelasticFS_F16" );
    // Class G4ParticleHPNT2AInelasticFS
    InsertModel( 25046, "model_G4ParticleHPNT2AInelasticFS_F17" );
    // Class G4ParticleHP4NInelasticFS
    InsertModel( 25047, "model_G4ParticleHP4NInelasticFS_F18" );
    // Class G4ParticleHP2NPInelasticFS
    InsertModel( 25048, "model_G4ParticleHP2NPInelasticFS_F19" );
    // Class G4ParticleHP3NPInelasticFS
    InsertModel( 25049, "model_G4ParticleHP3NPInelasticFS_F20" );
    // Class G4ParticleHPN2PInelasticFS
    InsertModel( 25050, "model_G4ParticleHPN2PInelasticFS_F21" );
    // Class G4ParticleHPNPAInelasticFS
    InsertModel( 25051, "model_G4ParticleHPNPAInelasticFS_F22" );
    // Class G4ParticleHPPInelasticFS
    InsertModel( 25052, "model_G4ParticleHPPInelasticFS_F23" );
    // Class G4ParticleHPDInelasticFS
    InsertModel( 25053, "model_G4ParticleHPDInelasticFS_F24" );
    // Class G4ParticleHPTInelasticFS
    InsertModel( 25054, "model_G4ParticleHPTInelasticFS_F25" );
    // Class G4ParticleHPHe3InelasticFS
    InsertModel( 25055, "model_G4ParticleHPHe3InelasticFS_F26" );
    // Class G4ParticleHPAInelasticFS
    InsertModel( 25056, "model_G4ParticleHPAInelasticFS_F27" );
    // Class G4ParticleHP2AInelasticFS
    InsertModel( 25057, "model_G4ParticleHP2AInelasticFS_F28" );
    // Class G4ParticleHP3AInelasticFS
    InsertModel( 25058, "model_G4ParticleHP3AInelasticFS_F29" );
    // Class G4ParticleHP2PInelasticFS
    InsertModel( 25059, "model_G4ParticleHP2PInelasticFS_F30" );
    // Class G4ParticleHPPAInelasticFS
    InsertModel( 25060, "model_G4ParticleHPPAInelasticFS_F31" );
    // Class G4ParticleHPD2AInelasticFS
    InsertModel( 25061, "model_G4ParticleHPD2AInelasticFS_F32" );
    // Class G4ParticleHPT2AInelasticFS
    InsertModel( 25062, "model_G4ParticleHPT2AInelasticFS_F33" );
    // Class G4ParticleHPPDInelasticFS
    InsertModel( 25063, "model_G4ParticleHPPDInelasticFS_F34" );
    // Class G4ParticleHPPTInelasticFS
    InsertModel( 25064, "model_G4ParticleHPPTInelasticFS_F35" );
    // Class G4ParticleHPDAInelasticFS
    InsertModel( 25065, "model_G4ParticleHPDAInelasticFS_F36" );

    // --- LEND : 25'200 - 25'299 ---

    // Class: G4LENDModel
    InsertModel( 25200, "model_LENDModel" );

    // The following classes inherit from G4LENDModel

    // Class: G4LENDCapture
    InsertModel( 25210, "model_LENDCapture" );
    // Class: G4LENDElastic
    InsertModel( 25220, "model_LENDElastic" );
    // Class: G4LENDFission
    InsertModel( 25230, "model_LENDFission" );
    // Class: G4LENDInelastic
    InsertModel( 25240, "model_LENDInelastic" );
    // Class: G4LENDCombinedModel
    InsertModel( 25250, "model_LENDCombinedModel" );
    // Class: G4LENDGammaModel
    InsertModel( 25260, "model_LENDGammaModel" );
    
    // --- Radioactive Decay : 25'500 - 25'999 ---
    
    // 25'510 +  10*G4RadioactiveDecayMode
    InsertModel( 25510, "model_RDM_IT" );
    InsertModel( 25520, "model_RDM_BetaMinus" );
    InsertModel( 25530, "model_RDM_BetaPlus" );
    InsertModel( 25540, "model_RDM_KshellEC" );
    InsertModel( 25550, "model_RDM_LshellEC" );
    InsertModel( 25560, "model_RDM_MshellEC" );
    InsertModel( 25570, "model_RDM_NshellEC" );
    InsertModel( 25580, "model_RDM_Alpha" );
    InsertModel( 25590, "model_RDM_Proton" );
    InsertModel( 25600, "model_RDM_Neutron" );
    InsertModel( 25610, "model_RDM_SpFission" );
    InsertModel( 25620, "model_RDM_BDProton" );
    InsertModel( 25630, "model_RDM_BDNeutron" );
    InsertModel( 25640, "model_RDM_Beta2Minus" );
    InsertModel( 25650, "model_RDM_Beta2Plus" );
    InsertModel( 25660, "model_RDM_Proton2" );
    InsertModel( 25670, "model_RDM_Neutron2" );
    InsertModel( 25680, "model_RDM_Triton" );

    InsertModel( 25810, "model_RDM_AtomicRelaxation" );

    // -------------------------------------------------------------------
    // --- Others HAD (everything not include above) : 26'000 - 26'999 ---
    // -------------------------------------------------------------------
    //    -  26'000 - 26'099 : Stopping
    //    -  26'100 - 26'199 : Fission
    //    -  26'200 - 26'299 : Abration
    //    -  26'300 - 26'399 : Coalescence
    //    ...

    // --- Stopping : 26'000 - 26'099 ---
    
    // Below are classes that derives from G4HadronStoppingProcess .
    // The names are the same as before Geant4 11, except for the "model_" prefix.
    // Classes that use it: G4HadronStoppingProcess .
    
    // Class: G4HadronicAbsorptionBertini
    InsertModel( 26000, "model_hBertiniCaptureAtRest_EMCascade" );
    InsertModel( 26001, "model_hBertiniCaptureAtRest_NuclearCapture" );
    InsertModel( 26002, "model_hBertiniCaptureAtRest_DIO" );
    // Class: G4HadronicAbsorptionFritiof
    InsertModel( 26010, "model_hFritiofCaptureAtRest_EMCascade" );
    InsertModel( 26011, "model_hFritiofCaptureAtRest_NuclearCapture" );
    InsertModel( 26012, "model_hFritiofCaptureAtRest_DIO" );
    // Class: G4HadronicAbsorptionFritiofWithBinaryCascade
    InsertModel( 26020, "model_hFritiofWithBinaryCascadeCaptureAtRest_EMCascade" );
    InsertModel( 26021, "model_hFritiofWithBinaryCascadeCaptureAtRest_NuclearCapture" );
    InsertModel( 26022, "model_hFritiofWithBinaryCascadeCaptureAtRest_DIO" );
    // Class: G4MuonMinusCapture
    InsertModel( 26030, "model_muMinusCaptureAtRest_EMCascade" );
    InsertModel( 26031, "model_muMinusCaptureAtRest_NuclearCapture" );
    InsertModel( 26032, "model_muMinusCaptureAtRest_DIO" );

    // --- Fission : 26'100 - 26'199 ---

    // Class G4LFission
    InsertModel( 26100, "model_G4LFission" );

    // LLNL fission (related classes: G4FissionLibrary, G4LLNLFission, G4FissLib, G4fissionEvent)
    InsertModel( 26110, "model_G4LLNLFission" );

    // --- Abration : 26'200 - 26'299 ---

    // Class G4WilsonAbrasionModel
    InsertModel( 26200, "model_G4WilsonAbrasion" );
    
    // --- Coalescence : 26'300 - 26'399 ---

    // Class G4CRCoalescence
    InsertModel( 26300, "model_G4CRCoalescence" );
    
    // ==========================================================================
    // === 3rd OTHER (i.e. non-EM and non-HAD) MODELS : from 30'000 to 39'999 ===
    // ==========================================================================

    // -------------------------------
    // --- Biasing : 30'000-30'999 ---
    // -------------------------------    

    // The name is the same as before Geant4 11, except for the "model_" prefix.
    // Classes that use it: G4BOptrForceCollision
    InsertModel( 30010, "model_GenBiasForceCollision" );

    // ----------------------------------
    // --- Channeling : 31'000-31'999 ---
    // ----------------------------------
    
    // The name is the same as before Geant4 11, except for the "model_" prefix.
    // Classes that use it: G4Channeling , G4ChannelingOptrChangeCrossSection
    InsertModel( 31010, "model_channeling" );

    // --- Others ... ---

    // ======================================================================
    // ================== 4th MODELS ADDED AFTER Geant4 11 ==================   
    // ======================================================================
    // PLEASE ADD MODELS ONLY BELOW HERE, WITH PROPER  modelID .
    // IF YOU ARE NOT SURE, PLEASE CONTACT ONE OF THE COORDINATORS OF THE
    // GEANT4 PHYSICS WORKING GROUPS.
    
    // ...
    
    SanityCheck();
    isInitialized = true;

    // The following call is commented out because it should be protected by
    // the verbosity level, but this would imply a dependency of the global
    // category to other categories which is not allowed.
    // Anyhow, the call  G4PhysicsModelCatalog::PrintAllInformation()
    // can be easily placed elsewhere (in either user code, or even in other
    // Geant4 classes).
    //PrintAllInformation();
  }
}

// --------------------------------------------------------------------------
void G4PhysicsModelCatalog::SanityCheck() {
  if ( theVectorOfModelIDs->size() != theVectorOfModelNames->size() ) {
    G4ExceptionDescription ed;
    ed << "theVectorOfModelIDs' size=" << theVectorOfModelIDs->size()
       << " is NOT the same as theVectorOfModelNames's size=" << theVectorOfModelNames->size();
    G4Exception( "G4PhysicsModelCatalog::SanityCheck()", "PhysModelCatalog001",
		 FatalException, ed, "should be the same!" );    
  } else {
    G4bool isModelIDOutsideRange = false;
    G4bool isModelIDRepeated = false;
    G4bool isModelNameRepeated = false;
    for ( int idx = 0; idx < Entries(); ++idx ) {
      G4int modelID = (*theVectorOfModelIDs)[ idx ];
      G4String modelName = (*theVectorOfModelNames)[ idx ];
      if ( modelID < GetMinAllowedModelIDValue() || modelID > GetMaxAllowedModelIDValue() ) {
	isModelIDOutsideRange = true;
      }
      for ( int jdx = idx + 1; jdx < Entries(); ++jdx ) {
        if(modelID == (*theVectorOfModelIDs)[jdx])
        {
          isModelIDRepeated = true;
        }
        if(modelName == (*theVectorOfModelNames)[jdx])
        {
          isModelNameRepeated = true;
        }
      }
    }
    if ( isModelIDOutsideRange || isModelIDRepeated || isModelNameRepeated ) {
      G4ExceptionDescription ed;
      if(isModelIDOutsideRange)
      {
        ed << "theVectorOfModelIDs has NOT all entries between "
           << GetMinAllowedModelIDValue() << " and "
           << GetMaxAllowedModelIDValue();
      }
      if(isModelIDRepeated)
      {
        ed << "theVectorOfModelIDs has NOT all unique IDs !";
      }
      if(isModelNameRepeated)
      {
        ed << "theVectorOfModelNames has NOT all unique names !";
      }
      G4Exception( "G4PhysicsModelCatalog::SanityCheck()", "PhysModelCatalog002",
		   FatalException, ed, "cannot continue!" );
    }
  }
  return;
}

// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
const G4String G4PhysicsModelCatalog::GetModelNameFromID( const G4int modelID ) {
  G4String modelName = "Undefined";
  if ( modelID >= GetMinAllowedModelIDValue()  &&  modelID <= GetMaxAllowedModelIDValue() ) {
    for ( int idx = 0; idx < Entries(); ++idx ) {
      if ( (*theVectorOfModelIDs)[ idx ] == modelID ) {
	modelName = (*theVectorOfModelNames)[ idx ];
	break;
      }
    }
  }
  return modelName;
}

// --------------------------------------------------------------------------
const G4String G4PhysicsModelCatalog::GetModelNameFromIndex( const G4int modelIndex ) {
  return ( modelIndex >= 0  &&  modelIndex < Entries() ) ? (*theVectorOfModelNames)[ modelIndex ] : "Undefined";
}

// --------------------------------------------------------------------------
G4int G4PhysicsModelCatalog::GetModelID( const G4int modelIndex ) {
  return ( modelIndex >= 0  &&  modelIndex < Entries() ) ? (*theVectorOfModelIDs)[ modelIndex ] : -1;
}

// --------------------------------------------------------------------------
G4int G4PhysicsModelCatalog::GetModelID( const G4String& modelName ) {
  if(!isInitialized)
  {
    Initialize();
  }
  G4int modelID = -1;
  for ( G4int idx = 0; idx < Entries(); ++idx ) {
    if ( (*theVectorOfModelNames)[ idx ] == modelName ) {
      modelID = (*theVectorOfModelIDs)[ idx ];
      break;
    }
  }
  return modelID;
}

// --------------------------------------------------------------------------
G4int G4PhysicsModelCatalog::GetModelIndex( const G4int modelID ) {
  G4int modelIndex = -1;
  if ( modelID >= GetMinAllowedModelIDValue()  &&  modelID <= GetMaxAllowedModelIDValue() ) {
    for ( G4int idx = 0; idx < Entries(); ++idx ) {
      if ( (*theVectorOfModelIDs)[ idx ] == modelID ) {
	modelIndex = idx;
	break;
      }
    }
  }
  return modelIndex;
}

// --------------------------------------------------------------------------
G4int G4PhysicsModelCatalog::GetModelIndex( const G4String& modelName ) {
  G4int modelIndex = -1;
  for ( G4int idx = 0; idx < Entries(); ++idx ) {
    if ( (*theVectorOfModelNames)[ idx ] == modelName ) {
      modelIndex = idx;
      break;
    }
  }
  return modelIndex;
}

// --------------------------------------------------------------------------
G4int G4PhysicsModelCatalog::Entries() {
  // It is enough to check the size of one of the two vectors, because they have the same size.
  return ( theVectorOfModelIDs != nullptr ) ? G4int( theVectorOfModelIDs->size() ) : -1;
}

// --------------------------------------------------------------------------
void G4PhysicsModelCatalog::PrintAllInformation() {
  G4cout << G4endl
	 << " ==================================================== " << G4endl
	 << " === G4PhysicsModelCatalog::PrintAllInformation() === " << G4endl
	 << " ==================================================== " << G4endl
	 << " SIZE (i.e. number of models in the catalog)=" << Entries() << G4endl;
  for ( int idx = 0; idx < Entries(); ++idx ) {
    G4int modelID = (*theVectorOfModelIDs)[ idx ];
    G4String modelName = (*theVectorOfModelNames)[ idx ];
    G4cout << "\t index=" << idx << "\t modelName=" << modelName
	   << "\t modelID=" << modelID << G4endl;
  }
  G4cout << " ==================================================== " << G4endl
	 << " ==================================================== " << G4endl
	 << " ==================================================== " << G4endl
	 << G4endl;
}
