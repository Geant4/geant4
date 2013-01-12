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
/// \file runAndEvent/RE02/include/RE02HadronPhysics.hh
/// \brief Definition of the RE02HadronPhysics class
//
// $Id$
// --------------------------------------------------------------
//
//  10-Oct-2003 Full Hadron Processes with Parameterization Model  T. Koi

#ifndef RE02HadronPhysics_h
#define RE02HadronPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"

// Hadronic Processes

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

// Binary Cascade Models
#include "G4BinaryCascade.hh"

// Low energy models
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

//
/// User hadron  physics constructor
///
///  applys related processes to hadrons
///
/// - void ConstructParticle()
///     does nothing
///
/// - void ConstructProcess()
///     adds processes to each particle
///     pi+ :
///       G4HadronElasticProcess with G4LElastic,
///       G4PionPlusInelasticProcess with G4LEPionPlusInelastic and
///       G4HEPionPlusInelastic,
///       G4hMultipleScattering and G4hIonisation
///     pi- :
///       G4HadronElasticProcess with G4LElastic,
///       G4PionMinusInelasticProcess with G4LEPionMinusInelastic and
///       G4HEPionMinusInelastic,
///       G4hMultipleScattering and G4hIonisation
///     K+ :
///       G4HadronElasticProcess with G4LElastic,
///       G4KaonPlusInelasticProcess with G4LEKaonPlusInelastic and
///       G4HEKaonPlusInelastic,
///       G4hMultipleScattering and G4hIonisation
///     K- :
///       G4HadronElasticProcess with G4LElastic,
///       G4KaonMinusInelasticProcess with G4LEKaonMinusInelastic and
///       G4HEKaonMinusInelastic,
///       G4hMultipleScattering and G4hIonisation
///     K0L :
///       G4HadronElasticProcess with G4LEastic and 
///       G4KaonZeoLInelasticProcess with G4LEKaonZeroLInelastic and
///       G4HEKaonZeroInelastic
///     K0S :
///       G4HadronElasticProcess with G4LEastic and
///       G4KaonZeroSInelasticProcess with G4LEKaonZeroSInelastic and
///       G4HEKaonZeroInelastic.
///     proton :
///       G4HadronElasticProcess with G4LElastic,
///       G4ProtonInelasticProcess with G4BinaryCascade, G4LEProtonInelastic and 
///       G4HEProtonInelastic,
///       G4hMultipleScattering and G4hIonisation
///     anti proton :
///       G4HadronElasticProcess with G4LElastic,
///       G4AntiProtonInelasticProcess with G4LEAntiProtonInelastic and 
///       G4HEAntiProtonInelastic,
///       G4AntiProtonAnnihilationAtRest, G4hMultipleScattering and G4hIonisation
///     neutron :
///       G4HadronElasticProcess with G4LElastic,
///       G4NeutronInelasticProcess with G4BinaryCascade, G4LENeutronInelastic
///       and G4HENeutronInelastic,
///       G4HadronFissionProcess with G4LFission and
///       G4HadronCaptureProcess with G4LCapture
///     anti neutron :
///       G4HadronElasticProcess with G4LElastic,
///       G4AntiNeutronInelasticProcess with G4LEAntiNeutronInelastic and 
///       G4HEAntiNeutronInelastic and
///       G4AntiNeutronAnnihilationAtRest
///     lambda :
///       G4HadronElasticProcess with G4LElastic,
///       G4LambdaInelasticProcess with G4LELambdaInelastic and 
///       G4HELambdaInelastic
///     anti lambda :
///       G4HadronElasticProcess with G4LElastic,
///       G4AntiLambdaInelasticProcess with G4LEAntiLambdaInelastic and 
///       G4HEAntiLambdaInelastic
///     sigma+ :
///       G4HadronElasticProcess with G4LElastic,
///       G4SigmaPlusInelasticProcess with G4LESigmaPlusInelastic and 
///       G4HESigmaPlusInelastic,
///       G4hMultipleScattering and G4hIonisation
///     anti sigma+ :
///       G4HadronElasticProcess with G4LElastic,
///       G4AntiSigmaPlusInelasticProcess with G4LEAntiSigmaPlusInelastic and 
///       G4HEAntiSigmaPlusInelastic,
///       G4hMultipleScattering and G4hIonisation
///     sigma- :
///       G4HadronElasticProcess with G4LElastic,
///       G4SigmaMinusInelasticProcess with G4LESigmaMinusInelastic and 
///       G4HESigmaMinusInelastic,
///       G4hMultipleScattering and G4hIonisation
///     anti sigma- :
///       G4HadronElasticProcess with G4LElastic,
///       G4AntiSigmaMinusInelasticProcess with G4LEAntiSigmaMinusInelastic and 
///       G4HEAntiSigmaMinusInelastic,
///       G4hMultipleScattering and G4hIonisation
///     Xi0 :
///       G4HadronElasticProcess with G4LElastic,
///       G4XiZeroInelasticProcess with G4LEXiZeroInelastic and 
///       G4HEXiZeroInelastic
///     anti Xi0 :
///       G4HadronElasticProcess with G4LElastic,
///       G4AntiXiZeroInelasticProcess with G4LEAntiXiZeroInelastic and 
///       G4HEAntiXiZeroInelastic
///     Xi- :
///       G4HadronElasticProcess with G4LElastic,
///       G4XiMinusInelasticProcess with G4LEXiMinusInelastic and 
///       G4HEXiMinusInelastic,
///       G4hMultipleScattering and G4hIonisation
///     anti Xi- :
///       G4HadronElasticProcess with G4LElastic,
///       G4AntiXiMinusInelasticProcess with G4LEAntiXiMinusInelastic and 
///       G4HEAntiXiMinusInelastic,
///       G4hMultipleScattering and G4hIonisation
///     omega- :
///       G4HadronElasticProcess with G4LElastic,
///       G4OmegaMinusInelasticProcess with G4LEOmegaMinusInelastic and 
///       G4HEOmegaMinusInelastic,
///       G4hMultipleScattering and G4hIonisation
///     anti omega- :
///       G4HadronElasticProcess with G4LElastic,
///       G4AntiOmegaMinusInelasticProcess with G4LEAntiOmegaMinusInelastic and 
///       G4HEAntiOmegaMinusInelastic,
///       G4hMultipleScattering and G4hIonisation
//
class RE02HadronPhysics : public G4VPhysicsConstructor
{
  public: 
    RE02HadronPhysics(const G4String& name="hadron");
    virtual ~RE02HadronPhysics();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle(){;};
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

  protected:

};

#endif

