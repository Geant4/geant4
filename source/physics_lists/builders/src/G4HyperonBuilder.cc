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
//---------------------------------------------------------------------------
// ClassName: G4HyperonBuilder
// Author: Alberto Ribon
// Date: May 2020
// Modified:
//---------------------------------------------------------------------------

#include "G4HyperonBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4HadronInelasticProcess.hh"


G4HyperonBuilder::G4HyperonBuilder() {  
  theLambdaInelastic         = new G4HadronInelasticProcess( "lambdaInelastic", G4Lambda::Definition() );
  theAntiLambdaInelastic     = new G4HadronInelasticProcess( "anti-lambdaInelastic", G4AntiLambda::Definition() );
  theSigmaMinusInelastic     = new G4HadronInelasticProcess( "sigma-Inelastic", G4SigmaMinus::Definition() );
  theAntiSigmaMinusInelastic = new G4HadronInelasticProcess( "anti_sigma-Inelastic", G4AntiSigmaMinus::Definition() );
  theSigmaPlusInelastic      = new G4HadronInelasticProcess( "sigma+Inelastic", G4SigmaPlus::Definition() );
  theAntiSigmaPlusInelastic  = new G4HadronInelasticProcess( "anti_sigma+Inelastic", G4AntiSigmaPlus::Definition() );
  theXiMinusInelastic        = new G4HadronInelasticProcess( "xi-Inelastic", G4XiMinus::Definition() );
  theAntiXiMinusInelastic    = new G4HadronInelasticProcess( "anti_xi-Inelastic", G4AntiXiMinus::Definition() );
  theXiZeroInelastic         = new G4HadronInelasticProcess( "xi0Inelastic", G4XiZero::Definition() );
  theAntiXiZeroInelastic     = new G4HadronInelasticProcess( "anti_xi0Inelastic", G4AntiXiZero::Definition() );
  theOmegaMinusInelastic     = new G4HadronInelasticProcess( "omega-Inelastic", G4OmegaMinus::Definition() );
  theAntiOmegaMinusInelastic = new G4HadronInelasticProcess( "anti_omega-Inelastic", G4AntiOmegaMinus::Definition() );
}


void G4HyperonBuilder::RegisterMe( G4PhysicsBuilderInterface* aB ) {
  auto bld = dynamic_cast< G4VHyperonBuilder* >( aB );
  if ( bld != nullptr ) theModelCollections.push_back( bld );
  else                  G4PhysicsBuilderInterface::RegisterMe( aB );
}


void G4HyperonBuilder::Build() {
  for ( std::vector< G4VHyperonBuilder* >::iterator i = theModelCollections.begin();
	i != theModelCollections.end(); ++i ) {
    (*i)->Build( theLambdaInelastic );
    (*i)->Build( theAntiLambdaInelastic );
    (*i)->Build( theSigmaMinusInelastic );
    (*i)->Build( theAntiSigmaMinusInelastic );
    (*i)->Build( theSigmaPlusInelastic );
    (*i)->Build( theAntiSigmaPlusInelastic );
    (*i)->Build( theXiMinusInelastic );
    (*i)->Build( theAntiXiMinusInelastic );
    (*i)->Build( theXiZeroInelastic );
    (*i)->Build( theAntiXiZeroInelastic );
    (*i)->Build( theOmegaMinusInelastic );
    (*i)->Build( theAntiOmegaMinusInelastic );
  }
  G4ProcessManager* aProcMan = nullptr;
  aProcMan = G4Lambda::Lambda()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theLambdaInelastic );
  aProcMan = G4AntiLambda::AntiLambda()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theAntiLambdaInelastic );
  aProcMan = G4SigmaMinus::SigmaMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theSigmaMinusInelastic );
  aProcMan = G4AntiSigmaMinus::AntiSigmaMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theAntiSigmaMinusInelastic );
  aProcMan = G4SigmaPlus::SigmaPlus()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theSigmaPlusInelastic );
  aProcMan = G4AntiSigmaPlus::AntiSigmaPlus()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theAntiSigmaPlusInelastic );
  aProcMan = G4XiMinus::XiMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theXiMinusInelastic );
  aProcMan = G4AntiXiMinus::AntiXiMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theAntiXiMinusInelastic );
  aProcMan = G4XiZero::XiZero()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theXiZeroInelastic );
  aProcMan = G4AntiXiZero::AntiXiZero()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theAntiXiZeroInelastic );
  aProcMan = G4OmegaMinus::OmegaMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theOmegaMinusInelastic );
  aProcMan = G4AntiOmegaMinus::AntiOmegaMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theAntiOmegaMinusInelastic );  
}
