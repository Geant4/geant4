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
// $Id: G4QNeutronBuilder.cc,v 1.4 2009-11-16 19:12:10 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4QNeutronBuilder
//
// Author: 2009 M. Kosov
//
// Modified:
//
//----------------------------------------------------------------------------
// Short description: for possible use in QGSC_QGSC (not needed for QGSC_CHIPS)
//-----------------------------------------------------------------------------
//
#include "G4QNeutronBuilder.hh"

G4QNeutronBuilder::G4QNeutronBuilder(): wasActivated(false) 
{
  theNeutronInelastic = new G4NeutronInelasticProcess;
  theCHIPSInelastic  = new G4QInelastic;
  const G4String& processName = "MixedNeutronInelasticProcess";
  const G4ParticleDefinition* proj = G4Neutron::Neutron();
  theProcessMixer= new G4QDiscProcessMixer(processName, proj);
}

G4QNeutronBuilder::~G4QNeutronBuilder() 
{
  delete theProcessMixer;
  delete theCHIPSInelastic;
  delete theNeutronInelastic;
}

void G4QNeutronBuilder::Build()
{
  wasActivated = true;
  std::vector<G4VNeutronBuilder *>::iterator i;
  for(i=theModelCollections.begin(); i!=theModelCollections.end(); i++)
  {
    (*i)->Build(theNeutronInelastic);
  }
  G4ProcessManager * theProcMan = G4Neutron::Neutron()->GetProcessManager();
  theProcessMixer->AddDiscreteProcess(theNeutronInelastic, 1.E8);// the second part is fake
  theProcessMixer->AddDiscreteProcess(theCHIPSInelastic, 290*megaelectronvolt);
  theProcMan->AddDiscreteProcess(theNeutronInelastic);
}
// 2009 by M. Kosov
