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
//      File name: G4ParticleHPInelasticDataPT.cc
//
//      Authors: Marek Zmeskal (CTU, Czech Technical University in Prague, Czech Republic)
//               Loic Thulliez (CEA France)      
//
//      Creation date: 4 June 2024
//
//      Description: Class for utilization of cross-sections from
//                   probability tables in the unresolved resonance region
//                   for inelastic channel.
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
//

#include "G4ParticleHPInelasticDataPT.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Neutron.hh"
#include "G4Element.hh"
#include "G4ParticleHPManager.hh"
#include "G4HadronicParameters.hh"
#include "G4ParticleHPProbabilityTablesStore.hh"
#include "G4HadronicException.hh"


G4ParticleHPInelasticDataPT::G4ParticleHPInelasticDataPT() : G4VCrossSectionDataSet( "NeutronHPInelasticXSPT" ) {
  // minimum and maximum energy limit for URR in ENDF/B-VII.1, it is overwritten in BuildPhysicsTable
  SetMinKinEnergy( 1.0 * CLHEP::eV );
  SetMaxKinEnergy( 1.2 * CLHEP::MeV );
  URRlimits = nullptr;
}


G4ParticleHPInelasticDataPT::~G4ParticleHPInelasticDataPT() {}


G4bool G4ParticleHPInelasticDataPT::IsIsoApplicable( const G4DynamicParticle* dp , G4int /*Z*/ , G4int /*A*/ ,
                                                     const G4Element* elm, const G4Material* /*mat*/ ) {
  // checks applicability for the element
  if ( doNOTusePTforInelastic ) {
      // do not use for njoy
      return false;
  }
  if ( dp->GetDefinition() != G4Neutron::Neutron() ) {
      return false;
    } else {
      std::size_t elementI = elm->GetIndex();
      G4double eKin = dp->GetKineticEnergy();
      if ( eKin < (*URRlimits).at(elementI).first ) {  // kinetic energy below the URR energy range for this element = minURR(isotopes in element)
        return false;
      } else if ( eKin > (*URRlimits).at(elementI).second ) {  // kinetic energy above the URR energy range for this element = maxURR(isotopes in element)
        return false;
      }
    return true;
  }
  return false;
}


G4double G4ParticleHPInelasticDataPT::GetIsoCrossSection( const G4DynamicParticle* dp, G4int /*Z*/ , G4int /*A*/ ,
                                                          const G4Isotope* iso, const G4Element* element, const G4Material* material ) {
  return G4ParticleHPProbabilityTablesStore::GetInstance()->GetIsoCrossSectionPT( dp, 3, iso, element, material );
}


void G4ParticleHPInelasticDataPT::BuildPhysicsTable( const G4ParticleDefinition& aP ) {
  if ( G4HadronicParameters::Instance()->GetTypeTablePT() == "njoy" ) {
    SetMinKinEnergy( DBL_MAX );
    SetMaxKinEnergy( 0.0 );
    doNOTusePTforInelastic = true;
  } else if ( G4HadronicParameters::Instance()->GetTypeTablePT() == "calendf" ) {
    doNOTusePTforInelastic = false;
    G4cout << "BuildPhysicsTable in G4ParticleHPInelasticDataPT." << G4endl;
    if ( &aP != G4Neutron::Neutron() ) {
        throw G4HadronicException( __FILE__, __LINE__, "Attempt to use NeutronHP data for particles other than neutrons!" );
    }
    URRlimits = G4ParticleHPManager::GetInstance()->GetURRlimits();
    if ( G4Threading::IsWorkerThread() ) {
      // sets the overall limits of the URR, which are stored at the last position of URRlimits - min and max URR(all elements)
      // defines URR model energy range
      SetMinKinEnergy( (*URRlimits).back().first );
      SetMaxKinEnergy( (*URRlimits).back().second );
    } else {
      if (G4ParticleHPManager::GetInstance()->GetProbabilityTables() == nullptr ) {
        G4ParticleHPProbabilityTablesStore::GetInstance()->Init();
        G4ParticleHPManager::GetInstance()->RegisterProbabilityTables( G4ParticleHPProbabilityTablesStore::GetInstance()->GetProbabilityTables() );
      }
      if ( URRlimits == nullptr ) {
        G4ParticleHPProbabilityTablesStore::GetInstance()->InitURRlimits();
        URRlimits = G4ParticleHPProbabilityTablesStore::GetInstance()->GetURRlimits();
        G4ParticleHPManager::GetInstance()->RegisterURRlimits( URRlimits );
      }
      // sets the overall limits of the URR, which are stored at the last position of URRlimits - min and max URR(all elements)
      // defines URR model energy range
      SetMinKinEnergy( (*URRlimits).back().first );
      SetMaxKinEnergy( (*URRlimits).back().second );
    }
  }
}


G4int G4ParticleHPInelasticDataPT::GetVerboseLevel() const {
  return G4ParticleHPManager::GetInstance()->GetVerboseLevel();
}


void G4ParticleHPInelasticDataPT::SetVerboseLevel( G4int newValue ) {
   G4ParticleHPManager::GetInstance()->SetVerboseLevel( newValue );
}


void G4ParticleHPInelasticDataPT::CrossSectionDescription( std::ostream& outFile ) const {
    outFile << "Inelastic probability tables." ;
}
