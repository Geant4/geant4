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
//
// ClassName:  G4HadronicAbsorptionFritiof
//
// Author:     Alberto Ribon
//
// Date:       27 July 2012
//
// Modified:
//
// Class Description:
//
// Intermediate class for hadronic absorption at rest using FTF/Preco.
// Physics lists should reference the concrete subclasses for:
// anti_proton, anti_sigma+, and all anti-nuclei.
//
//---------------------------------------------------------------------------

#include <iostream>

#include "G4SystemOfUnits.hh"
#include "G4HadronicAbsorptionFritiof.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4TheoFSGenerator.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4HadronicInteractionRegistry.hh"

// Constructor
G4HadronicAbsorptionFritiof::
G4HadronicAbsorptionFritiof( G4ParticleDefinition* pdef )
  : G4HadronStoppingProcess( "hFritiofCaptureAtRest" ), 
    pdefApplicable( pdef ) {
  
  G4TheoFSGenerator * theModel = new G4TheoFSGenerator( "FTFP" );
  G4FTFModel * theStringModel = new G4FTFModel;
  theLund = new G4LundStringFragmentation;
  theStringDecay = new G4ExcitedStringDecay( theLund );
  theStringModel->SetFragmentationModel( theStringDecay );

  // Not a cascade - goes straight to Preco
  G4HadronicInteraction* p =
    G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
  G4VPreCompoundModel * thePreEquilib = static_cast<G4VPreCompoundModel*>(p); 
  if(! thePreEquilib) { thePreEquilib = new G4PreCompoundModel; }

  G4GeneratorPrecompoundInterface * theCascade = 
    new G4GeneratorPrecompoundInterface( thePreEquilib ); 

  theModel->SetHighEnergyGenerator( theStringModel );
  theModel->SetTransport( theCascade );

  G4double theMin = 0.0*GeV;
  G4double theMax = 100.0*TeV;
  theModel->SetMinEnergy( theMin );
  theModel->SetMaxEnergy( theMax );

  RegisterMe( theModel );
}


G4HadronicAbsorptionFritiof::~G4HadronicAbsorptionFritiof() {
  delete theLund;
  delete theStringDecay;
}


// Applies to constructor-specified particle, or to all known cases
G4bool G4HadronicAbsorptionFritiof::
IsApplicable( const G4ParticleDefinition& particle ) {
  return ( ( 0 == pdefApplicable && 
             ( &particle == G4AntiProton::Definition() ||
	       &particle == G4AntiSigmaPlus::Definition() ||
               particle.GetBaryonNumber() < -1 ) )     // Anti-nuclei
	   || ( &particle == pdefApplicable ) );
}


// Documentation of purpose
void G4HadronicAbsorptionFritiof::
ProcessDescription( std::ostream& os ) const {
  os << "Stopping and absorption of anti_protons, anti_sigma+, and \n"
     << "all anti-nuclei using Fritiof (FTF) string model.\n"
     << "Geant4 PreCompound model is used for nuclear de-excitation."
     << std::endl;
}
