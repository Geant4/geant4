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
//    **************************************
//    *                                    *
//    *    HadrontherapyProtonHadro.hh        *
//    *                                    *
//    **************************************
//
// $Id: HadrontherapyProtonHadro.hh,v 1.5 2005-05-27 15:05:53 cirrone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author : Susanna Guatelli, guatelli@ge.infn.it
// 

#ifndef HadrontherapyProtonHadro_h
#define HadrontherapyProtonHadro_h 1

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

class HadrontherapyProtonHadro: public G4VPhysicsConstructor 
{
 public:
   HadrontherapyProtonHadro(const G4String& name = "proton-hadronic");
   virtual ~HadrontherapyProtonHadro();

 protected:
   // Construct particle and physics
   void ConstructParticle(){};
   void ConstructProcess();

 private:
  G4ProtonInelasticCrossSection theXSec;
  G4ExcitationHandler theHandler;  
  
  G4DeuteronInelasticProcess      theIPdeuteron;
  G4TritonInelasticProcess        theIPtriton;
  G4AlphaInelasticProcess         theIPalpha;
  G4HadronInelasticProcess*       theIPHe3;
  G4HadronInelasticProcess*       theIPIonC12;
  G4HadronInelasticProcess*       theIPGenericIon;
  G4ProtonInelasticCrossSection  thePXSec;
  G4ProtonInelasticProcess       theIPProton;
  G4NeutronInelasticProcess	 theIPNeutron;
  G4NeutronInelasticCrossSection   theNXSec;
};
#endif








