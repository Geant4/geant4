// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BarionConstructor.cc,v 1.1 1999-01-07 16:10:00 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//

#include "G4BarionConstructor.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
// Barions
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Neutron.hh"
#include "G4AntiNeutron.hh"

#include "G4Lambda.hh"
#include "G4SigmaPlus.hh"
#include "G4SigmaZero.hh"
#include "G4SigmaMinus.hh"
#include "G4XiMinus.hh"
#include "G4XiZero.hh"
#include "G4OmegaMinus.hh"

#include "G4AntiLambda.hh"
#include "G4AntiSigmaPlus.hh"
#include "G4AntiSigmaZero.hh"
#include "G4AntiSigmaMinus.hh"
#include "G4AntiXiMinus.hh"
#include "G4AntiXiZero.hh"
#include "G4AntiOmegaMinus.hh"

#include "G4LambdacPlus.hh"
#include "G4SigmacPlusPlus.hh"
#include "G4SigmacPlus.hh"
#include "G4SigmacZero.hh"
#include "G4XicPlus.hh"
#include "G4XicZero.hh"
#include "G4OmegacZero.hh"

#include "G4AntiLambdacPlus.hh"
#include "G4AntiSigmacPlusPlus.hh"
#include "G4AntiSigmacPlus.hh"
#include "G4AntiSigmacZero.hh"
#include "G4AntiXicPlus.hh"
#include "G4AntiXicZero.hh"
#include "G4AntiOmegacZero.hh"

G4BarionConstructor::G4BarionConstructor()
{
}

G4BarionConstructor::~G4BarionConstructor()
{
}


void G4BarionConstructor::ConstructParticle()
{
  ConstructNucleons();
  ConstructStrangeBarions();
  ConstructCharmBarions();
  ConstructBottomBarions();
}

void G4BarionConstructor::ConstructNucleons()
{
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
}
void G4BarionConstructor::ConstructStrangeBarions()
{
  G4Lambda::LambdaDefinition();
  G4AntiLambda::AntiLambdaDefinition();
  G4SigmaZero::SigmaZeroDefinition();
  G4AntiSigmaZero::AntiSigmaZeroDefinition();
  G4SigmaPlus::SigmaPlusDefinition();
  G4AntiSigmaPlus::AntiSigmaPlusDefinition();
  G4SigmaMinus::SigmaMinusDefinition();
  G4AntiSigmaMinus::AntiSigmaMinusDefinition();
  G4XiZero::XiZeroDefinition();
  G4AntiXiZero::AntiXiZeroDefinition();
  G4XiMinus::XiMinusDefinition();
  G4AntiXiMinus::AntiXiMinusDefinition();
  G4OmegaMinus::OmegaMinusDefinition();
  G4AntiOmegaMinus::AntiOmegaMinusDefinition();
}
void G4BarionConstructor::ConstructCharmBarions()
{
  G4LambdacPlus::LambdacPlusDefinition();
  G4SigmacPlusPlus::SigmacPlusPlusDefinition();
  G4SigmacPlus::SigmacPlusDefinition();
  G4SigmacZero::SigmacZeroDefinition();
  G4XicPlus::XicPlusDefinition();
  G4XicZero::XicZeroDefinition();
  G4OmegacZero::OmegacZeroDefinition();
  G4AntiLambdacPlus::AntiLambdacPlusDefinition();
  G4AntiSigmacPlusPlus::AntiSigmacPlusPlusDefinition();
  G4AntiSigmacPlus::AntiSigmacPlusDefinition();
  G4AntiSigmacZero::AntiSigmacZeroDefinition();
  G4AntiXicPlus::AntiXicPlusDefinition();
  G4AntiXicZero::AntiXicZeroDefinition();
  G4AntiOmegacZero::AntiOmegacZeroDefinition();
}

void G4BarionConstructor::ConstructBottomBarions()
{
}
