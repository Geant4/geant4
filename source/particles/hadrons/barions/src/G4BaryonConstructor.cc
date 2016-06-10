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
// $Id: G4BaryonConstructor.cc 67971 2013-03-13 10:13:24Z gcosmo $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//

#include "G4BaryonConstructor.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
// Baryons
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

#include "G4Lambdab.hh"
#include "G4SigmabPlus.hh"
#include "G4SigmabZero.hh"
#include "G4SigmabMinus.hh"
#include "G4XibZero.hh"
#include "G4XibMinus.hh"
#include "G4OmegabMinus.hh"

#include "G4AntiLambdab.hh"
#include "G4AntiSigmabPlus.hh"
#include "G4AntiSigmabZero.hh"
#include "G4AntiSigmabMinus.hh"
#include "G4AntiXibZero.hh"
#include "G4AntiXibMinus.hh"
#include "G4AntiOmegabMinus.hh"


G4BaryonConstructor::G4BaryonConstructor()
{
}

G4BaryonConstructor::~G4BaryonConstructor()
{
}


void G4BaryonConstructor::ConstructParticle()
{
  ConstructNucleons();
  ConstructStrangeBaryons();
  ConstructCharmBaryons();
  ConstructBottomBaryons();
}

void G4BaryonConstructor::ConstructNucleons()
{
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
}
void G4BaryonConstructor::ConstructStrangeBaryons()
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
void G4BaryonConstructor::ConstructCharmBaryons()
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

void G4BaryonConstructor::ConstructBottomBaryons()
{
  G4Lambdab::LambdabDefinition();
  G4SigmabPlus::SigmabPlusDefinition();
  G4SigmabZero::SigmabZeroDefinition();
  G4SigmabMinus::SigmabMinusDefinition();
  G4XibZero::XibZeroDefinition();
  G4XibMinus::XibMinusDefinition();
  G4OmegabMinus::OmegabMinusDefinition();

  G4AntiLambdab::AntiLambdabDefinition();
  G4AntiSigmabPlus::AntiSigmabPlusDefinition();
  G4AntiSigmabZero::AntiSigmabZeroDefinition();
  G4AntiSigmabMinus::AntiSigmabMinusDefinition();
  G4AntiXibZero::AntiXibZeroDefinition();
  G4AntiXibMinus::AntiXibMinusDefinition();
  G4AntiOmegabMinus::AntiOmegabMinusDefinition();

}
