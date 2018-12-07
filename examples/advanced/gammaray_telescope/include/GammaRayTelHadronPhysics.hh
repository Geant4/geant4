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
#include "G4NeutronCaptureXS.hh"
#include "G4HadronElastic.hh"

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

// Low-energy Models
#include "G4LFission.hh"

// Cross section handlers and high energy models

#include "G4VCrossSectionDataSet.hh"
#include "G4CascadeInterface.hh"
#include "G4CrossSectionPairGG.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4ChipsKaonZeroInelasticXS.hh"
#include "G4ChipsKaonMinusInelasticXS.hh"
#include "G4ChipsKaonPlusInelasticXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4BGGNucleonInelasticXS.hh"

// Stopping processes
#include "G4AntiProtonAbsorptionFritiof.hh"
#include "G4PiMinusAbsorptionBertini.hh"
#include "G4KaonMinusAbsorptionBertini.hh"

// quark gluon string model with chips afterburner.
#include "G4FTFModel.hh"
#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
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
  G4HadronElastic*        theElasticModel;
  
   // Pi + 
   G4PionPlusInelasticProcess thePionPlusInelastic;
   G4hMultipleScattering thePionPlusMult;
   G4hIonisation thePionPlusIonisation;

   // Pi -
   G4PionMinusInelasticProcess thePionMinusInelastic;
   G4hMultipleScattering thePionMinusMult;
   G4hIonisation thePionMinusIonisation;
   G4PiMinusAbsorptionBertini thePionMinusAbsorption;

   // pi+ and pi-
   
    G4TheoFSGenerator* theModel;
    G4ExcitationHandler theHandler;
    G4PreCompoundModel * thePreEquilib;
    G4GeneratorPrecompoundInterface* theCascade;
    G4FTFModel* theStringModel;
    G4QGSMFragmentation theFragmentation;
    G4ExcitedStringDecay * theStringDecay;

   // K + 
   G4KaonPlusInelasticProcess theKaonPlusInelastic;
   G4hMultipleScattering theKaonPlusMult;
   G4hIonisation theKaonPlusIonisation;
	
   // K -
   G4KaonMinusInelasticProcess theKaonMinusInelastic;
   G4hMultipleScattering theKaonMinusMult;
   G4hIonisation theKaonMinusIonisation;
   G4KaonMinusAbsorptionBertini theKaonMinusAbsorption;

   // K0L
   G4KaonZeroLInelasticProcess theKaonZeroLInelastic;

   // K0S
   G4KaonZeroSInelasticProcess theKaonZeroSInelastic;

   // Proton
   G4ProtonInelasticProcess theProtonInelastic; 
   G4hMultipleScattering theProtonMult;
   G4hIonisation theProtonIonisation;
 
   // anti-proton
   G4AntiProtonInelasticProcess theAntiProtonInelastic; 
   G4hMultipleScattering theAntiProtonMult;
   G4hIonisation theAntiProtonIonisation;
   G4AntiProtonAbsorptionFritiof  theAntiProtonAnnihilation;
    
   // neutron
   G4NeutronInelasticProcess  theNeutronInelastic; 
   G4HadronFissionProcess theNeutronFission;
   G4LFission* theNeutronFissionModel;
   G4HadronCaptureProcess*  theNeutronCapture;

   // anti-neutron
   G4AntiNeutronInelasticProcess  theAntiNeutronInelastic;  
   
  
};


#endif

