//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// ********************************************************************

#ifndef LISAPhysicsList_h
#define LISAPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

// Hadronics
#include "G4HadronElasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4HadronFissionProcess.hh"

// Inelastic Processes
#include "G4ProtonInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4LambdaInelasticProcess.hh"
#include "G4AntiLambdaInelasticProcess.hh"
#include "G4OmegaMinusInelasticProcess.hh"
#include "G4AntiOmegaMinusInelasticProcess.hh"
#include "G4SigmaMinusInelasticProcess.hh"
#include "G4AntiSigmaMinusInelasticProcess.hh"
#include "G4SigmaPlusInelasticProcess.hh"
#include "G4AntiSigmaPlusInelasticProcess.hh"
#include "G4XiZeroInelasticProcess.hh"
#include "G4AntiXiZeroInelasticProcess.hh"
#include "G4XiMinusInelasticProcess.hh"
#include "G4AntiXiMinusInelasticProcess.hh"

// Low-energy Parameterised Models: 1 to 25 GeV
#include "G4LElastic.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4LEKaonZeroSInelastic.hh"
#include "G4LEKaonZeroLInelastic.hh"
#include "G4LEKaonMinusInelastic.hh"
// #include "G4LEProtonInelastic.hh"
#include "G4LEAntiProtonInelastic.hh"
// #include "G4LENeutronInelastic.hh"
#include "G4LEAntiNeutronInelastic.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEAlphaInelastic.hh"
#include "G4LELambdaInelastic.hh"
#include "G4LEAntiLambdaInelastic.hh"
#include "G4LEOmegaMinusInelastic.hh"
#include "G4LEAntiOmegaMinusInelastic.hh"
#include "G4LESigmaMinusInelastic.hh"
#include "G4LEAntiSigmaMinusInelastic.hh"
#include "G4LESigmaPlusInelastic.hh"
#include "G4LEAntiSigmaPlusInelastic.hh"
#include "G4LEXiZeroInelastic.hh"
#include "G4LEAntiXiZeroInelastic.hh"
#include "G4LEXiMinusInelastic.hh"
#include "G4LEAntiXiMinusInelastic.hh"
// neutrons
#include "G4LCapture.hh"
#include "G4LFission.hh"

// High-energy Parameterised Models: 25 GeV to 10 TeV
//  #include "G4HEPionPlusInelastic.hh"
//  #include "G4HEPionMinusInelastic.hh"
//  #include "G4HEKaonPlusInelastic.hh"
//  #include "G4HEKaonZeroInelastic.hh"
//  #include "G4HEKaonZeroInelastic.hh"
//  #include "G4HEKaonMinusInelastic.hh"
//  #include "G4HEProtonInelastic.hh"
//  #include "G4HENeutronInelastic.hh"
#include "G4HEAntiProtonInelastic.hh"
#include "G4HEAntiNeutronInelastic.hh"
#include "G4HELambdaInelastic.hh"
#include "G4HEAntiLambdaInelastic.hh"
#include "G4HEOmegaMinusInelastic.hh"
#include "G4HEAntiOmegaMinusInelastic.hh"
#include "G4HESigmaMinusInelastic.hh"
#include "G4HEAntiSigmaMinusInelastic.hh"
#include "G4HESigmaPlusInelastic.hh"
#include "G4HEAntiSigmaPlusInelastic.hh"
#include "G4HEXiZeroInelastic.hh"
#include "G4HEAntiXiZeroInelastic.hh"
#include "G4HEXiMinusInelastic.hh"
#include "G4HEAntiXiMinusInelastic.hh"

// // Neutron HP Models: Thermal to 19 MeV
// #include "G4NeutronHPElastic.hh"
// #include "G4NeutronHPElasticData.hh"
// #include "G4NeutronHPCapture.hh"
// #include "G4NeutronHPCaptureData.hh"
// #include "G4NeutronHPInelastic.hh"
// #include "G4NeutronHPInelasticData.hh"

// Stopping processes
#include "G4PiMinusAbsorptionAtRest.hh"
#include "G4KaonMinusAbsorptionAtRest.hh"
#include "G4AntiProtonAnnihilationAtRest.hh"
#include "G4AntiNeutronAnnihilationAtRest.hh"

// Generator models: HE
#include "G4TheoFSGenerator.hh"
#include "G4Evaporation.hh"
#include "G4CompetitiveFission.hh"
#include "G4FermiBreakUp.hh"
#include "G4StatMF.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

// Kinetic Model
#include "G4BinaryCascade.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"


///////////////////////////
// ElectroNuclear Physics

// photonuclear and electronuclear reaction
#include "G4PhotoNuclearProcess.hh"
#include "G4ElectronNuclearProcess.hh"
#include "G4PositronNuclearProcess.hh"
#include "G4GammaNuclearReaction.hh"
#include "G4ElectroNuclearReaction.hh"

// CHIPS fragmentation model
#include "G4TheoFSGenerator.hh"
#include "G4StringChipsParticleLevelInterface.hh"
#include "G4QGSModel.hh"
#include "G4GammaParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

// muon photonuclear reaction
#include "G4MuNuclearInteraction.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class LISAPhysicsList: public G4VUserPhysicsList {

  public:
    LISAPhysicsList();
    virtual ~LISAPhysicsList();

  public:
    virtual void SetCuts();

  protected:
    // particles and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();
    
    // physics processes
    virtual void AddTransportation();
    virtual void ElectromagneticPhysics();
    virtual void HadronicPhysics();
    virtual void ElectroNuclearPhysics();
    virtual void GeneralPhysics();

  private:
    G4int VerboseLevel;
  protected:

  // PhotoNuclear *************************************************

  G4PhotoNuclearProcess thePhotoNuclearProcess;
  G4GammaNuclearReaction* theGammaReaction;      
  G4TheoFSGenerator* theHEModel_PN;
  G4StringChipsParticleLevelInterface* theCascade_PN;
  G4QGSModel<G4GammaParticipants> theStringModel_PN;
  G4QGSMFragmentation theFragmentation_PN;
  G4ExcitedStringDecay* theStringDecay_PN;
  // ElectronNuclear
  G4ElectronNuclearProcess theElectronNuclearProcess;
  G4ElectroNuclearReaction* theElectroReaction;
  // PositronNuclear
  G4PositronNuclearProcess thePositronNuclearProcess;
  // MuNucleus
  G4MuNuclearInteraction theMuMinusNuclearInteraction;
  G4MuNuclearInteraction theMuPlusNuclearInteraction;
  

   // Hadronics  *************************************************
  
  // Binary Cascade
  G4TheoFSGenerator* theHEModel;
  G4Evaporation* theEvaporation;
  G4FermiBreakUp* theFermiBreakUp;
  G4StatMF* theMF;
  G4ExcitationHandler* theHandler;
  G4PreCompoundModel* thePreEquilib;
  G4GeneratorPrecompoundInterface* theCascade;
  G4VPartonStringModel* theStringModel;
  G4BinaryCascade* theCasc;
  G4VLongitudinalStringDecay* theFragmentation;
  G4ExcitedStringDecay* theStringDecay;
  G4BinaryCascade* theCascForPi;
  // Cascade for light ions
  G4BinaryLightIonReaction* theIonCascade;
  G4TripathiCrossSection* theTripathiCrossSection;
  G4IonsShenCrossSection* theShenCrossSection;
  G4BinaryLightIonReaction* theGenIonCascade;
  
  
   // Elastic Process
  G4HadronElasticProcess theElasticProcess;
  G4LElastic* theElasticModel;

  // pi+
  G4PionPlusInelasticProcess thePionPlusInelasticProcess;
  G4LEPionPlusInelastic* theLEPionPlusInelasticModel;

  // pi-
  G4PionMinusInelasticProcess thePionMinusInelasticProcess;
  G4LEPionMinusInelastic* theLEPionMinusInelasticModel;
  G4PiMinusAbsorptionAtRest thePiMinusAbsorptionAtRest;
  
  // kaon+
  G4KaonPlusInelasticProcess theKaonPlusInelasticProcess;
  G4LEKaonPlusInelastic* theLEKaonPlusInelasticModel;

  // kaon0S
  G4KaonZeroSInelasticProcess theKaonZeroSInelasticProcess;
  G4LEKaonZeroSInelastic* theLEKaonZeroSInelasticModel;

  // kaon0L
  G4KaonZeroLInelasticProcess theKaonZeroLInelasticProcess;
  G4LEKaonZeroLInelastic* theLEKaonZeroLInelasticModel;

  // kaon-
  G4KaonMinusInelasticProcess theKaonMinusInelasticProcess;
  G4LEKaonMinusInelastic* theLEKaonMinusInelasticModel;
  G4KaonMinusAbsorptionAtRest theKaonMinusAbsorptionAtRest;
  
  // proton
  G4ProtonInelasticProcess theProtonInelasticProcess;

  // anti-proton
  G4AntiProtonInelasticProcess theAntiProtonInelasticProcess;
  G4LEAntiProtonInelastic* theLEAntiProtonInelasticModel;
  G4HEAntiProtonInelastic* theHEAntiProtonInelasticModel;
  G4AntiProtonAnnihilationAtRest theAntiProtonAnnihilationAtRest;

  // neutron
  G4HadronElasticProcess theNeutronElasticProcess;
  G4LElastic* theNeutronElasticModel1;
  //   G4NeutronHPElastic* theNeutronElasticModel2;
  //   G4NeutronHPElasticData* theNeutronElasticData;
  G4NeutronInelasticProcess theNeutronInelasticProcess;
  //   G4NeutronHPInelastic* theNeutronInelasticModel1;
  //   G4NeutronHPInelasticData* theNeutronInelasticData1;
  G4HadronCaptureProcess theNeutronCaptureProcess;
  G4LCapture* theNeutronCaptureModel1;
  //   G4NeutronHPCapture* theNeutronCaptureModel2;
  //   G4NeutronHPCaptureData* theNeutronCaptureData;
  G4HadronFissionProcess theNeutronFissionProcess;
  G4LFission* theNeutronFissionModel;

  // anti-neutron
  G4AntiNeutronInelasticProcess theAntiNeutronInelasticProcess;
  G4LEAntiNeutronInelastic* theLEAntiNeutronInelasticModel;
  G4HEAntiNeutronInelastic* theHEAntiNeutronInelasticModel;
  G4AntiNeutronAnnihilationAtRest theAntiNeutronAnnihilationAtRest;

  // deuteron
  G4DeuteronInelasticProcess* theDeuteronInelasticProcess;
  G4LEDeuteronInelastic* theLEDeuteronInelasticModel;

  // triton
  G4TritonInelasticProcess* theTritonInelasticProcess;
  G4LETritonInelastic* theLETritonInelasticModel;

  // alpha
  G4AlphaInelasticProcess* theAlphaInelasticProcess;
  G4LEAlphaInelastic* theLEAlphaInelasticModel;

  // He-3
  G4HadronInelasticProcess* theHe3InelasticProcess;

  // Generic Ion
  G4HadronInelasticProcess* theGenericIonInelasticProcess;

  // lambda
  G4LambdaInelasticProcess theLambdaInelasticProcess;
  G4LELambdaInelastic* theLELambdaInelasticModel;
  G4HELambdaInelastic* theHELambdaInelasticModel;

  // anti-lambda
  G4AntiLambdaInelasticProcess theAntiLambdaInelasticProcess;
  G4LEAntiLambdaInelastic* theLEAntiLambdaInelasticModel;
  G4HEAntiLambdaInelastic* theHEAntiLambdaInelasticModel;

  // omega-
  G4OmegaMinusInelasticProcess theOmegaMinusInelasticProcess;
  G4LEOmegaMinusInelastic* theLEOmegaMinusInelasticModel;
  G4HEOmegaMinusInelastic* theHEOmegaMinusInelasticModel;

  // anti-omega-
  G4AntiOmegaMinusInelasticProcess theAntiOmegaMinusInelasticProcess;
  G4LEAntiOmegaMinusInelastic* theLEAntiOmegaMinusInelasticModel;
  G4HEAntiOmegaMinusInelastic* theHEAntiOmegaMinusInelasticModel;

  // sigma-
  G4SigmaMinusInelasticProcess theSigmaMinusInelasticProcess;
  G4LESigmaMinusInelastic* theLESigmaMinusInelasticModel;
  G4HESigmaMinusInelastic* theHESigmaMinusInelasticModel;

  // anti-sigma-
  G4AntiSigmaMinusInelasticProcess theAntiSigmaMinusInelasticProcess;
  G4LEAntiSigmaMinusInelastic* theLEAntiSigmaMinusInelasticModel;
  G4HEAntiSigmaMinusInelastic* theHEAntiSigmaMinusInelasticModel;

  // sigma+
  G4SigmaPlusInelasticProcess theSigmaPlusInelasticProcess;
  G4LESigmaPlusInelastic* theLESigmaPlusInelasticModel;
  G4HESigmaPlusInelastic* theHESigmaPlusInelasticModel;

  // anti-sigma+
  G4AntiSigmaPlusInelasticProcess theAntiSigmaPlusInelasticProcess;
  G4LEAntiSigmaPlusInelastic* theLEAntiSigmaPlusInelasticModel;
  G4HEAntiSigmaPlusInelastic* theHEAntiSigmaPlusInelasticModel;

  // xi0
  G4XiZeroInelasticProcess theXiZeroInelasticProcess;
  G4LEXiZeroInelastic* theLEXiZeroInelasticModel;
  G4HEXiZeroInelastic* theHEXiZeroInelasticModel;

  // anti-xi0
  G4AntiXiZeroInelasticProcess theAntiXiZeroInelasticProcess;
  G4LEAntiXiZeroInelastic* theLEAntiXiZeroInelasticModel;
  G4HEAntiXiZeroInelastic* theHEAntiXiZeroInelasticModel;

  // xi-
  G4XiMinusInelasticProcess theXiMinusInelasticProcess;
  G4LEXiMinusInelastic* theLEXiMinusInelasticModel;
  G4HEXiMinusInelastic* theHEXiMinusInelasticModel;

  // anti-xi-
  G4AntiXiMinusInelasticProcess theAntiXiMinusInelasticProcess;
  G4LEAntiXiMinusInelastic* theLEAntiXiMinusInelasticModel;
  G4HEAntiXiMinusInelastic* theHEAntiXiMinusInelasticModel;


};

#endif
