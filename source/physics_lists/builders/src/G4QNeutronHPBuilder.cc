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
//---------------------------------------------------------------------------
//
// ClassName:   G4QNeutronHPBuilder
//
// Author: April 2012 M. Kosov
//
// Modified:
//
//----------------------------------------------------------------------------
// Short description: for use in CHIPS_HP physics list (mix Q w/ HP at 20 MeV)
//----------------------------------------------------------------------------
//
#include "G4QNeutronHPBuilder.hh"
#include "G4SystemOfUnits.hh"

G4QNeutronHPBuilder::G4QNeutronHPBuilder(): 
    theNeutrons(0)
    , theNeutronFission(0)
    , theNeutronCapture(0)
    , theCHIPSNGamma(0)
    , theHPNeutron(0)
    , wasActivated(false) 
{
  theNeutronInelastic = new G4NeutronInelasticProcess;
  theCHIPSInelastic  = new G4QInelastic;
  const G4ParticleDefinition* proj = G4Neutron::Neutron();
  const G4String& INprocessName = "MixedHPQNeutronInelasticProcess";
  theInProcessMixer= new G4QDiscProcessMixer(INprocessName, proj);
  const G4String& NGprocessName = "MixedHPQNeutronGammaProcess";
  theNgProcessMixer= new G4QDiscProcessMixer(NGprocessName, proj);
  const G4String& FIprocessName = "MixedHPQNeutronFissionProcess";
  theFiProcessMixer= new G4QDiscProcessMixer(FIprocessName, proj);
}

G4QNeutronHPBuilder::~G4QNeutronHPBuilder() 
{
  delete theCHIPSInelastic;
  delete theCHIPSNGamma;
  delete theNeutronInelastic;
  delete theHPNeutron;
  delete theInProcessMixer;
  delete theNgProcessMixer;
  delete theFiProcessMixer;
}

void G4QNeutronHPBuilder::Build()
{
  wasActivated = true;
  // The model definition for neutrons (needed for HP implementation)
  theNeutrons = new G4NeutronBuilder;
  theNeutrons->RegisterMe( theHPNeutron = new G4NeutronHPBuilder );
  // End of the model definition
  std::vector<G4VNeutronBuilder *>::iterator i;
  for(i = theModelCollections.begin(); i != theModelCollections.end(); i++)
  {
    (*i)->Build(theNeutronInelastic);
    (*i)->Build(theNeutronCapture);
    (*i)->Build(theNeutronFission);
  }
  G4ProcessManager * theProcMan = G4Neutron::Neutron()->GetProcessManager();

  theInProcessMixer->AddDiscreteProcess(theCHIPSInelastic, 1.E8*megaelectronvolt);
  theInProcessMixer->AddDiscreteProcess(theNeutronInelastic, 19.9*megaelectronvolt);

  theNgProcessMixer->AddDiscreteProcess(theCHIPSNGamma, 1.E8*megaelectronvolt);
  theNgProcessMixer->AddDiscreteProcess(theNeutronCapture, 19.9*megaelectronvolt);

  //theFiProcessMixer->AddDiscreteProcess(theCHIPSFission, 1.E8*megaelectronvolt);
  //theFiProcessMixer->AddDiscreteProcess(theNeutronFission, 19.9*megaelectronvolt);

  theProcMan->AddDiscreteProcess(theInProcessMixer); // Mix CHIPS+HP for inelastic
  theProcMan->AddDiscreteProcess(theNgProcessMixer); // Mix CHIPS+HP for (n,gamma)
  theProcMan->AddDiscreteProcess(theNeutronFission); // Only HP for fission
  //theProcMan->AddDiscreteProcess(theFiProcessMixer); // Mix CHIPS+HP for fission
}
// 2009 by M. Kosov
