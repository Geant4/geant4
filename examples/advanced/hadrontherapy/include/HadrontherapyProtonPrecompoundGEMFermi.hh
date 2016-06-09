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
// $Id: HadrontherapyProtonPrecompoundGEMFermi.hh; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------
#ifndef HadrontherapyProtonPrecompoundGEMFermi_h
#define HadrontherapyProtonPrecompoundGEMFermi_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4ProtonInelasticProcess.hh"

#include "G4ProcessManager.hh"

#include "G4ExcitationHandler.hh"  
#include "G4PreCompoundModel.hh"   

#include "G4ProtonInelasticCrossSection.hh"
#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4ExcitedStringDecay.hh"
// For neutron
#include "G4HENeutronInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"
// For ions
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4BinaryCascade.hh"
#include "G4IonInelasticProcess.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4LETritonInelastic.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4LEAlphaInelastic.hh"

class HadrontherapyProtonPrecompoundGEMFermi: public G4VPhysicsConstructor 
{
 public:
   HadrontherapyProtonPrecompoundGEMFermi(const G4String& name = "proton-precompoundGEMFermi");
   virtual ~HadrontherapyProtonPrecompoundGEMFermi();

 protected:
   // Construct particle and physics
   void ConstructParticle(){};
   void ConstructProcess();

 private:

  G4double binaryLightIonLowLimit;
  G4double binaryLightIonHighLimit;
  G4double LEPHighLimit;
  G4double precompoundLowLimit;
  G4double precompoundHighLimit;
  G4double neutronLowLimit;
  G4double neutronHighLimit;
  G4int targetZ;
  G4int targetA;

  G4ExcitationHandler theHandler; 


  // Proton inelastic proces
  G4ProtonInelasticProcess       theIPProton;
  // Cross Section for proton inelastic process
  G4ProtonInelasticCrossSection  thePXSec;

  // Neutron inelastic process
  G4NeutronInelasticProcess	 theIPNeutron;
  // Cross Section for neutron inelastic process
  G4NeutronInelasticCrossSection   theNXSec;

  // Deuteron inelastic process
  G4DeuteronInelasticProcess      theIPdeuteron;
     
  // Tritium inelastic process
  G4TritonInelasticProcess        theIPtriton;

  // Alpha inelastic process
  G4AlphaInelasticProcess         theIPalpha;
 
  // He3 inelastic process
  G4HadronInelasticProcess*       theIPHe3;
};
#endif








