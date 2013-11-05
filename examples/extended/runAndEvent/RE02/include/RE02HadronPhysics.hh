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

// Bertini Cascade model
#include "G4CascadeInterface.hh"

// Low energy models
#include "G4HadronElastic.hh"
#include "G4LFission.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronCaptureXS.hh"

// High energy models
#include "G4TheoFSGenerator.hh"
#include "G4FTFModel.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4LundStringFragmentation.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4PreCompoundModel.hh"
  
// Stopping processes
#include "G4AntiProtonAbsorptionFritiof.hh"
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
///       G4HadronElasticProcess with G4HadronElastic,
///       G4PionPlusInelasticProcess with G4CascadeInterface (Bertini) and FTF,
///       G4hMultipleScattering and G4hIonisation
///     pi- :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4PionMinusInelasticProcess with G4CascadeInterface (Bertini) and FTF,
///       G4hMultipleScattering and G4hIonisation
///     K+ :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4KaonPlusInelasticProcess with G4CascadeInterface (Bertini) and FTF,
///       G4hMultipleScattering and G4hIonisation
///     K- :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4KaonMinusInelasticProcess with G4CascadeInterface and FTF,
///       G4hMultipleScattering and G4hIonisation
///     K0L :
///       G4HadronElasticProcess with G4HadronElastic and 
///       G4KaonZeoLInelasticProcess with G4CascadeInterface and FTF
///     K0S :
///       G4HadronElasticProcess with G4HadronElastic and
///       G4KaonZeroSInelasticProcess with G4CascadeInterface and FTF
///     proton :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4ProtonInelasticProcess with G4BinaryCascade, G4CascadeInterface and 
///       FTF, 
///       G4hMultipleScattering and G4hIonisation
///     anti proton :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4AntiProtonInelasticProcess with FTF,
///       G4AntiProtonAbsorptionFritiof, G4hMultipleScattering and G4hIonisation
///     neutron :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4NeutronInelasticProcess with G4BinaryCascade, G4CascadeInterface and
///       FTF,
///       G4HadronFissionProcess with G4LFission and
///       G4HadronCaptureProcess with G4NeutronRadCapture
///     anti neutron :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4AntiNeutronInelasticProcess with FTF,
///       G4AntiNeutronAnnihilationAtRest
///     lambda :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4LambdaInelasticProcess with G4CascadeInterface and FTF
///     anti lambda :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4AntiLambdaInelasticProcess with FTF
///     sigma+ :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4SigmaPlusInelasticProcess with G4LCascadeInterface and FTF,
///       G4hMultipleScattering and G4hIonisation
///     anti sigma+ :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4AntiSigmaPlusInelasticProcess with FTF, 
///       G4hMultipleScattering and G4hIonisation
///     sigma- :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4SigmaMinusInelasticProcess with G4CascadeInterface and FTF, 
///       G4hMultipleScattering and G4hIonisation
///     anti sigma- :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4AntiSigmaMinusInelasticProcess with FTF
///       G4hMultipleScattering and G4hIonisation
///     Xi0 :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4XiZeroInelasticProcess with G4CascadeInterface and FTF
///     anti Xi0 :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4AntiXiZeroInelasticProcess with FTF
///     Xi- :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4XiMinusInelasticProcess with G4CascadeInterface and FTF,
///       G4hMultipleScattering and G4hIonisation
///     anti Xi- :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4AntiXiMinusInelasticProcess with FTF
///       G4hMultipleScattering and G4hIonisation
///     omega- :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4OmegaMinusInelasticProcess with G4CascadeInterface and FTF
///       G4hMultipleScattering and G4hIonisation
///     anti omega- :
///       G4HadronElasticProcess with G4HadronElastic,
///       G4AntiOmegaMinusInelasticProcess with FTF
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

