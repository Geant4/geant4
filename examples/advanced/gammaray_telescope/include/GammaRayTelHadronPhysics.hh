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
// $Id$
//
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This class is an derived class of G4VPhysicsConstructor
//
// ------------------------------------------------------------
#ifndef GammaRayTelHadronPhysics_h
#define GammaRayTelHadronPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <vector>

#include "G4VPhysicsConstructor.hh"

#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"

#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4LambdaInelasticProcess.hh"
#include "G4AntiLambdaInelasticProcess.hh"
#include "G4SigmaPlusInelasticProcess.hh"
#include "G4SigmaMinusInelasticProcess.hh"
#include "G4AntiSigmaPlusInelasticProcess.hh"
#include "G4AntiSigmaMinusInelasticProcess.hh"
#include "G4XiZeroInelasticProcess.hh"
#include "G4XiMinusInelasticProcess.hh"
#include "G4AntiXiZeroInelasticProcess.hh"
#include "G4AntiXiMinusInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4OmegaMinusInelasticProcess.hh"
#include "G4AntiOmegaMinusInelasticProcess.hh"

// Low-energy Models
#include "G4LElastic.hh"   
#include "G4LFission.hh"
#include "G4LCapture.hh"

#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4LEKaonZeroSInelastic.hh"
#include "G4LEKaonZeroLInelastic.hh"
#include "G4LEKaonMinusInelastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LEAntiProtonInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4LEAntiNeutronInelastic.hh"
#include "G4LELambdaInelastic.hh"
#include "G4LEAntiLambdaInelastic.hh"
#include "G4LESigmaPlusInelastic.hh"
#include "G4LESigmaMinusInelastic.hh"
#include "G4LEAntiSigmaPlusInelastic.hh"
#include "G4LEAntiSigmaMinusInelastic.hh"
#include "G4LEXiZeroInelastic.hh"
#include "G4LEXiMinusInelastic.hh"
#include "G4LEAntiXiZeroInelastic.hh"
#include "G4LEAntiXiMinusInelastic.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEAlphaInelastic.hh"
#include "G4LEOmegaMinusInelastic.hh"
#include "G4LEAntiOmegaMinusInelastic.hh"

// High-energy Models

#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4HEKaonPlusInelastic.hh"
#include "G4HEKaonZeroInelastic.hh"
#include "G4HEKaonZeroInelastic.hh"
#include "G4HEKaonMinusInelastic.hh"
#include "G4HEProtonInelastic.hh"
#include "G4HEAntiProtonInelastic.hh"
#include "G4HENeutronInelastic.hh"
#include "G4HEAntiNeutronInelastic.hh"
#include "G4HELambdaInelastic.hh"
#include "G4HEAntiLambdaInelastic.hh"
#include "G4HESigmaPlusInelastic.hh"
#include "G4HESigmaMinusInelastic.hh"
#include "G4HEAntiSigmaPlusInelastic.hh"
#include "G4HEAntiSigmaMinusInelastic.hh"
#include "G4HEXiZeroInelastic.hh"
#include "G4HEXiMinusInelastic.hh"
#include "G4HEAntiXiZeroInelastic.hh"
#include "G4HEAntiXiMinusInelastic.hh"
#include "G4HEOmegaMinusInelastic.hh"
#include "G4HEAntiOmegaMinusInelastic.hh"

// Stopping processes
#include "G4AntiProtonAnnihilationAtRest.hh"
#include "G4AntiNeutronAnnihilationAtRest.hh"

#ifdef TRIUMF_STOP_PIMINUS
#include "G4PionMinusAbsorptionAtRest.hh"
#else
#include "G4PiMinusAbsorptionAtRest.hh"
#endif
#ifdef TRIUMF_STOP_KMINUS
#include "G4KaonMinusAbsorption.hh"
#else
#include "G4KaonMinusAbsorptionAtRest.hh"
#endif

// quark gluon string model with chips afterburner.
#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"


class GammaRayTelHadronPhysics : public G4VPhysicsConstructor
{
  public: 
    GammaRayTelHadronPhysics(const G4String& name ="hadron");
    virtual ~GammaRayTelHadronPhysics();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    void ConstructParticle(){};
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
  void ConstructProcess();

  protected:
   // Elastic Process
    G4HadronElasticProcess theElasticProcess;
    G4LElastic*            theElasticModel;
  
   // Pi + 
   G4PionPlusInelasticProcess thePionPlusInelastic;
   G4LEPionPlusInelastic* theLEPionPlusModel;
   G4hMultipleScattering thePionPlusMult;
   G4hIonisation thePionPlusIonisation;

   // Pi -
   G4PionMinusInelasticProcess thePionMinusInelastic;
   G4LEPionMinusInelastic* theLEPionMinusModel;
   G4hMultipleScattering thePionMinusMult;
   G4hIonisation thePionMinusIonisation;
#ifdef TRIUMF_STOP_PIMINUS
   G4PionMinusAbsorptionAtRest thePionMinusAbsorption;
#else
   G4PiMinusAbsorptionAtRest thePionMinusAbsorption;
#endif

   // pi+ and pi-
   
    G4TheoFSGenerator* theModel;
    G4ExcitationHandler theHandler;
    G4PreCompoundModel * thePreEquilib;
    G4GeneratorPrecompoundInterface* theCascade;
    G4QGSModel<G4QGSParticipants>* theStringModel;
    G4QGSMFragmentation theFragmentation;
    G4ExcitedStringDecay * theStringDecay;

   // K + 
   G4KaonPlusInelasticProcess theKaonPlusInelastic;
   G4LEKaonPlusInelastic* theLEKaonPlusModel;
   G4HEKaonPlusInelastic* theHEKaonPlusModel;
   G4hMultipleScattering theKaonPlusMult;
   G4hIonisation theKaonPlusIonisation;
	
   // K -
   G4KaonMinusInelasticProcess theKaonMinusInelastic;
   G4LEKaonMinusInelastic* theLEKaonMinusModel;
   G4HEKaonMinusInelastic* theHEKaonMinusModel;
   G4hMultipleScattering theKaonMinusMult;
   G4hIonisation theKaonMinusIonisation;
#ifdef TRIUMF_STOP_KMINUS
   G4KaonMinusAbsorption theKaonMinusAbsorption;
#else
   G4PiMinusAbsorptionAtRest theKaonMinusAbsorption;
#endif

   // K0L
   G4KaonZeroLInelasticProcess theKaonZeroLInelastic;
   G4LEKaonZeroLInelastic* theLEKaonZeroLModel;
   G4HEKaonZeroInelastic* theHEKaonZeroLModel;

   // K0S
   G4KaonZeroSInelasticProcess theKaonZeroSInelastic;
   G4LEKaonZeroSInelastic* theLEKaonZeroSModel;
   G4HEKaonZeroInelastic* theHEKaonZeroSModel;

   // Proton
   G4ProtonInelasticProcess theProtonInelastic;
   G4LEProtonInelastic* theLEProtonModel;
   G4HEProtonInelastic* theHEProtonModel;
   G4hMultipleScattering theProtonMult;
   G4hIonisation theProtonIonisation;
 
   // anti-proton
   G4AntiProtonInelasticProcess theAntiProtonInelastic;
   G4LEAntiProtonInelastic* theLEAntiProtonModel;
   G4HEAntiProtonInelastic* theHEAntiProtonModel;
   G4hMultipleScattering theAntiProtonMult;
   G4hIonisation theAntiProtonIonisation;
   G4AntiProtonAnnihilationAtRest  theAntiProtonAnnihilation;
    
   // neutron
   G4NeutronInelasticProcess  theNeutronInelastic;
   G4LENeutronInelastic* theLENeutronModel;
   G4HENeutronInelastic* theHENeutronModel;
   G4HadronFissionProcess theNeutronFission;
   G4LFission* theNeutronFissionModel;
   G4HadronCaptureProcess  theNeutronCapture;
   G4LCapture* theNeutronCaptureModel;


   // anti-neutron
   G4AntiNeutronInelasticProcess  theAntiNeutronInelastic;
   G4LEAntiNeutronInelastic* theLEAntiNeutronModel;
   G4HEAntiNeutronInelastic* theHEAntiNeutronModel;
   G4AntiNeutronAnnihilationAtRest  theAntiNeutronAnnihilation;
   
   // Lambda
   G4LambdaInelasticProcess  theLambdaInelastic;
   G4LELambdaInelastic*  theLELambdaModel;
   G4HELambdaInelastic*  theHELambdaModel;
  
   // AntiLambda
   G4AntiLambdaInelasticProcess  theAntiLambdaInelastic;
   G4LEAntiLambdaInelastic*  theLEAntiLambdaModel;
   G4HEAntiLambdaInelastic*  theHEAntiLambdaModel;
  
   // SigmaMinus
   G4SigmaMinusInelasticProcess  theSigmaMinusInelastic;
   G4LESigmaMinusInelastic*  theLESigmaMinusModel;
   G4HESigmaMinusInelastic*  theHESigmaMinusModel;
   G4hMultipleScattering theSigmaMinusMult;
   G4hIonisation theSigmaMinusIonisation;
  
   // AntiSigmaMinus
   G4AntiSigmaMinusInelasticProcess  theAntiSigmaMinusInelastic;
   G4LEAntiSigmaMinusInelastic*  theLEAntiSigmaMinusModel;
   G4HEAntiSigmaMinusInelastic*  theHEAntiSigmaMinusModel;
   G4hMultipleScattering theAntiSigmaMinusMult;
   G4hIonisation theAntiSigmaMinusIonisation;
   
   // SigmaPlus
   G4SigmaPlusInelasticProcess  theSigmaPlusInelastic;
   G4LESigmaPlusInelastic*  theLESigmaPlusModel;
   G4HESigmaPlusInelastic*  theHESigmaPlusModel;
   G4hMultipleScattering theSigmaPlusMult;
   G4hIonisation theSigmaPlusIonisation;
  
   // AntiSigmaPlus
   G4AntiSigmaPlusInelasticProcess  theAntiSigmaPlusInelastic;
   G4LEAntiSigmaPlusInelastic*  theLEAntiSigmaPlusModel;
   G4HEAntiSigmaPlusInelastic*  theHEAntiSigmaPlusModel;
   G4hMultipleScattering theAntiSigmaPlusMult;
   G4hIonisation theAntiSigmaPlusIonisation;
  
   // XiZero
   G4XiZeroInelasticProcess  theXiZeroInelastic;
   G4LEXiZeroInelastic*  theLEXiZeroModel;
   G4HEXiZeroInelastic*  theHEXiZeroModel;
  
   // AntiXiZero
   G4AntiXiZeroInelasticProcess  theAntiXiZeroInelastic;
   G4LEAntiXiZeroInelastic*  theLEAntiXiZeroModel;
   G4HEAntiXiZeroInelastic*  theHEAntiXiZeroModel;
  
   // XiMinus
   G4XiMinusInelasticProcess  theXiMinusInelastic;
   G4LEXiMinusInelastic*  theLEXiMinusModel;
   G4HEXiMinusInelastic*  theHEXiMinusModel;
   G4hMultipleScattering theXiMinusMult;
   G4hIonisation theXiMinusIonisation;

   // AntiXiMinus
   G4AntiXiMinusInelasticProcess  theAntiXiMinusInelastic;
   G4LEAntiXiMinusInelastic*  theLEAntiXiMinusModel;
   G4HEAntiXiMinusInelastic*  theHEAntiXiMinusModel;
   G4hMultipleScattering theAntiXiMinusMult;
   G4hIonisation theAntiXiMinusIonisation;
  
   // OmegaMinus
   G4OmegaMinusInelasticProcess  theOmegaMinusInelastic;
   G4LEOmegaMinusInelastic*  theLEOmegaMinusModel;
   G4HEOmegaMinusInelastic*  theHEOmegaMinusModel;
   G4hMultipleScattering theOmegaMinusMult;
   G4hIonisation theOmegaMinusIonisation;
   
   // AntiOmegaMinus
   G4AntiOmegaMinusInelasticProcess  theAntiOmegaMinusInelastic;
   G4LEAntiOmegaMinusInelastic*  theLEAntiOmegaMinusModel;
   G4HEAntiOmegaMinusInelastic*  theHEAntiOmegaMinusModel;
   G4hMultipleScattering theAntiOmegaMinusMult;
   G4hIonisation theAntiOmegaMinusIonisation;   
};


#endif

