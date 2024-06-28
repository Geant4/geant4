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
// -------------------------------------------------------------------
//
//      Geant4 source file 
//
//      File name: G4ParticleHPCaptureURR.cc
//
//      Authors: Marek Zmeskal (CTU, Czech Technical University in Prague, Czech Republic)
//               Loic Thulliez (CEA France)      
//
//      Creation date: 4 June 2024
//
//      Description: Class to handle URR range, can be omitted once the
//                   proper isotope cross-section is stored in ParticleHP.
//
//      Modifications:
//      
// -------------------------------------------------------------------
// 
//

#include "G4ParticleHPCaptureURR.hh"
#include "G4ParticleHPManager.hh"
#include "G4ParticleHPChannel.hh"
#include "G4NeutronHPCapture.hh"
#include "G4ParticleHPProbabilityTablesStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"


G4ParticleHPCaptureURR::G4ParticleHPCaptureURR() : G4HadronicInteraction( "NeutronHPCaptureURR" ) {
  SetMinEnergy(  0.0 * CLHEP::eV );
  SetMaxEnergy( 20.0 * CLHEP::MeV );
  neutronHPcapture = new G4NeutronHPCapture;
}


G4ParticleHPCaptureURR::~G4ParticleHPCaptureURR() {}


G4HadFinalState* G4ParticleHPCaptureURR::ApplyYourself( const G4HadProjectile& aTrack, G4Nucleus& aNucleus ) {
  const G4Material* theMaterial = aTrack.GetMaterial();
  G4double kineticEnergy = aTrack.GetKineticEnergy();
  G4HadFinalState* theFinalState = nullptr;
  if ( kineticEnergy < (*URRlimits).back().first  ||  kineticEnergy >(*URRlimits).back().second ) {
    return neutronHPcapture->ApplyYourself( aTrack, aNucleus );
  }
  G4int elementI = -1;
  G4int isotopeJ = -1;
  G4int A = aNucleus.GetA_asInt();
  G4int Z = aNucleus.GetZ_asInt();
  // finds the element and isotope of the selected target aNucleus
  for ( G4int i = 0; i < (G4int)theMaterial->GetNumberOfElements(); ++i ) {
    if ( Z == theMaterial->GetElement(i)->GetZasInt() ) {
      for ( G4int j = 0; j < (G4int)theMaterial->GetElement(i)->GetNumberOfIsotopes(); ++j ) {
        if ( A == theMaterial->GetElement(i)->GetIsotope(j)->GetN() ) {
          isotopeJ = j;
          break;
        }
      }
      // the loop cannot be ended here because the material can have two elements with same Z but different isotopic composition
      if ( isotopeJ != -1 ) {
        // isotope was found and for loop is ended
        elementI = (G4int)theMaterial->GetElement(i)->GetIndex();
        break;
      }
    }  // end if find element
  }  // end element loop
  // Check whether the energy is out of the URR limits for the given element
  if ( kineticEnergy < (*URRlimits).at(elementI).first  ||  kineticEnergy > (*URRlimits).at(elementI).second ) { 
    // Call capture final state in G4ParicleHPChannel and SELECT ISOTOPE (to be improved in the future)
    G4ParticleHPManager::GetInstance()->OpenReactionWhiteBoard();
    theFinalState = (*G4ParticleHPManager::GetInstance()->GetCaptureFinalStates())[elementI]->ApplyYourself( aTrack );
    // Update target nucleus information according to the selected isotope
    G4int selectedIsotope_A = G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargA();
    aNucleus.SetParameters( selectedIsotope_A, Z );
    const G4Element* target_element = (*G4Element::GetElementTable())[elementI];
    const G4Isotope* target_isotope = nullptr;
    // Find the selected isotope among in the element
    for ( G4int j = 0; j < (G4int)target_element->GetNumberOfIsotopes(); ++j ) {
      target_isotope = target_element->GetIsotope(j);
      if ( target_isotope->GetN() == selectedIsotope_A ) break;
    }
    aNucleus.SetIsotope( target_isotope );
    G4ParticleHPManager::GetInstance()->CloseReactionWhiteBoard();
  } else {  // the energy is inside the limits of the URR
    // calls the final state for the found element and isotope
    theFinalState = ((*G4ParticleHPManager::GetInstance()->GetCaptureFinalStates())[elementI]->GetFinalStates())[isotopeJ]->ApplyYourself( aTrack );
  }
  return theFinalState;
}


void G4ParticleHPCaptureURR::BuildPhysicsTable( const G4ParticleDefinition& ) {
  neutronHPcapture->BuildPhysicsTable( *(G4Neutron::Neutron() ) );
  URRlimits = G4ParticleHPManager::GetInstance()->GetURRlimits();
  if ( URRlimits == nullptr ) {
    G4ParticleHPProbabilityTablesStore::GetInstance()->InitURRlimits();
    URRlimits = G4ParticleHPProbabilityTablesStore::GetInstance()->GetURRlimits();
    G4ParticleHPManager::GetInstance()->RegisterURRlimits( URRlimits );
  }
}


const std::pair< G4double, G4double > G4ParticleHPCaptureURR::GetFatalEnergyCheckLevels() const {
  // max energy non-conservation is mass of heavy nucleus
  return std::pair< G4double, G4double >( 10.0 * perCent, 350.0 * CLHEP::GeV );
}


G4int G4ParticleHPCaptureURR::GetVerboseLevel() const {
  return G4ParticleHPManager::GetInstance()->GetVerboseLevel();
}


void G4ParticleHPCaptureURR::SetVerboseLevel( G4int newValue ) {
  G4ParticleHPManager::GetInstance()->SetVerboseLevel( newValue );
}


void G4ParticleHPCaptureURR::ModelDescription( std::ostream& outFile ) const {
   outFile << "High Precision model based on Evaluated Nuclear Data Files (ENDF) for capture reaction of neutrons in the unresolved resonance region.";
}
