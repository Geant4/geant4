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
// --------------------------------------------------------------
//	GEANT 4 class implementation file
//

#include "G4IonConstructor.hh"
// Nuclei
#include "G4Alpha.hh"
#include "G4Deuteron.hh"
#include "G4GenericIon.hh"
#include "G4He3.hh"
#include "G4Triton.hh"
// AntiNuclei
#include "G4AntiAlpha.hh"
#include "G4AntiDeuteron.hh"
#include "G4AntiHe3.hh"
#include "G4AntiTriton.hh"
// Hyper-nuclei
#include "G4DoubleHyperDoubleNeutron.hh"
#include "G4DoubleHyperH4.hh"
#include "G4HyperAlpha.hh"
#include "G4HyperH4.hh"
#include "G4HyperHe5.hh"
#include "G4HyperTriton.hh"
// Anti-hyper-nuclei
#include "G4AntiDoubleHyperDoubleNeutron.hh"
#include "G4AntiDoubleHyperH4.hh"
#include "G4AntiHyperAlpha.hh"
#include "G4AntiHyperH4.hh"
#include "G4AntiHyperHe5.hh"
#include "G4AntiHyperTriton.hh"

void G4IonConstructor::ConstructParticle()
{
  ConstructLightIons();
  ConstructHyperNuclei();
}

void G4IonConstructor::ConstructLightIons()
{
  //  nuclei
  G4Alpha::AlphaDefinition();
  G4Deuteron::DeuteronDefinition();
  G4Triton::TritonDefinition();
  G4He3::He3Definition();
  //  anti_nuclei
  G4AntiAlpha::AntiAlphaDefinition();
  G4AntiDeuteron::AntiDeuteronDefinition();
  G4AntiTriton::AntiTritonDefinition();
  G4AntiHe3::AntiHe3Definition();
  //  generic ion
  G4GenericIon::GenericIonDefinition();
}

void G4IonConstructor::ConstructHyperNuclei()
{
  G4DoubleHyperDoubleNeutron::DoubleHyperDoubleNeutron();
  G4DoubleHyperH4::DoubleHyperH4();
  G4HyperAlpha::HyperAlpha();
  G4HyperH4::HyperH4();
  G4HyperHe5::HyperHe5();
  G4HyperTriton::HyperTriton();

  G4AntiDoubleHyperDoubleNeutron::AntiDoubleHyperDoubleNeutron();
  G4AntiDoubleHyperH4::AntiDoubleHyperH4();
  G4AntiHyperAlpha::AntiHyperAlpha();
  G4AntiHyperH4::AntiHyperH4();
  G4AntiHyperHe5::AntiHyperHe5();
  G4AntiHyperTriton::AntiHyperTriton();
}
