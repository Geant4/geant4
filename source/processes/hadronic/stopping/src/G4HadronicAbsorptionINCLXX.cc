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
// ClassName:  G4HadronicAbsorptionINCLXX
//
// Author:     Alberto Ribon
//
// Date:       April 2023
//
// Modified:
//
// Class Description:
//
// Intermediate class for hadronic absorption at rest using INCLXX/Preco.
//
//---------------------------------------------------------------------------

#include <iostream>

#include "G4SystemOfUnits.hh"
#include "G4HadronicAbsorptionINCLXX.hh"
#include "G4PreCompoundModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4HadronicParameters.hh"
#include "G4INCLXXInterface.hh" 


G4HadronicAbsorptionINCLXX::G4HadronicAbsorptionINCLXX( G4ParticleDefinition* pdef ) : 
  G4HadronStoppingProcess( "hINCLXXCaptureAtRest" ), pdefApplicable( pdef ) 
{
  G4INCLXXInterface* theModel = new G4INCLXXInterface;
  G4double theMin =  0.0*GeV;
  G4double theMax = 10.0*GeV;
  theModel->SetMinEnergy( theMin );
  theModel->SetMaxEnergy( theMax );
  RegisterMe( theModel );
}


G4HadronicAbsorptionINCLXX::~G4HadronicAbsorptionINCLXX() {}


// Applies to constructor-specified particle, or to all known cases
G4bool G4HadronicAbsorptionINCLXX::IsApplicable( const G4ParticleDefinition& particle ) {
  return ( (  pdefApplicable == nullptr  && ( &particle == G4AntiProton::Definition() ) )
	   || ( &particle == pdefApplicable ) );
}


// Documentation of purpose
void G4HadronicAbsorptionINCLXX::ProcessDescription( std::ostream& os ) const {
  os << "Stopping and absorption of anti_proton using INCLXX model.\n"
     << "Geant4 PreCompound model is used for nuclear de-excitation."
     << std::endl;
}
