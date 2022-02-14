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
// ClassName:  G4HadronicAbsorptionFritiofWithBinaryCascade
//
// Author:     Alberto Ribon
//
// Date:       July 2019
//
// Modified:
//
// Class Description:
//
// Intermediate class for hadronic absorption at rest using FTF/BIC.
//
//---------------------------------------------------------------------------

#include <iostream>

#include "G4SystemOfUnits.hh"
#include "G4HadronicAbsorptionFritiofWithBinaryCascade.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4TheoFSGenerator.hh"
#include "G4BinaryCascade.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4HadronicParameters.hh"

// Constructor
G4HadronicAbsorptionFritiofWithBinaryCascade::
G4HadronicAbsorptionFritiofWithBinaryCascade( G4ParticleDefinition* pdef )
  : G4HadronStoppingProcess( "hFritiofWithBinaryCascadeCaptureAtRest" ), 
    pdefApplicable( pdef ) {
  
  G4TheoFSGenerator* theModel = new G4TheoFSGenerator( "FTFB" );
  G4FTFModel* theStringModel = new G4FTFModel;
  theStringModel->SetFragmentationModel(new G4ExcitedStringDecay());

  G4BinaryCascade* theCascade = new G4BinaryCascade;

  theModel->SetHighEnergyGenerator( theStringModel );
  theModel->SetTransport( theCascade );

  G4double theMin = 0.0*GeV;
  G4double theMax = G4HadronicParameters::Instance()->GetMaxEnergy();
  theModel->SetMinEnergy( theMin );
  theModel->SetMaxEnergy( theMax );

  RegisterMe( theModel );
}

G4HadronicAbsorptionFritiofWithBinaryCascade::~G4HadronicAbsorptionFritiofWithBinaryCascade() {
}

// Applies to constructor-specified particle, or to all known cases
G4bool G4HadronicAbsorptionFritiofWithBinaryCascade::
IsApplicable( const G4ParticleDefinition& particle ) {
  // It is applicable only for anti_proton and anti_neutron, but not for
  // light anti-ion because the propagate method for nucleus-nucleus interaction:
  //   G4VIntraNuclearTransportModel::Propagate()
  // is not implemented.
  return ( ( ( ! pdefApplicable )  && 
             ( &particle == G4AntiProton::Definition()  ||  &particle == G4AntiNeutron::Definition() ) )
	   ||  &particle == pdefApplicable );
}


// Documentation of purpose
void G4HadronicAbsorptionFritiofWithBinaryCascade::
ProcessDescription( std::ostream& os ) const {
  os << "Stopping and absorption of anti_proton and anti_neutron \n"
     << "using  Fritiof (FTF) string model coupled with Binary Cascade (BIC) model."
     << std::endl;
}
