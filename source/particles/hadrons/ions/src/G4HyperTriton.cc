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
// Author:       2021 Alberto Ribon
//
//----------------------------------------------------------------------------

#include "G4HyperTriton.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                      HYPERTRITON                               ###
// ######################################################################

G4HyperTriton* G4HyperTriton::theInstance = nullptr;


G4HyperTriton* G4HyperTriton::Definition() {
  if ( theInstance != nullptr ) return theInstance;
  const G4String name = "hypertriton";
  // search in particle table
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4Ions* anInstance =  static_cast< G4Ions* >( pTable->FindParticle( name ) );
  if ( anInstance == nullptr ) {
    // create particle
    //
    //    Arguments for constructor are as follows
    //                            name             mass           width         charge
    //                          2*spin           parity   C-conjugation
    //                       2*Isospin       2*Isospin3        G-parity
    //                            type    lepton number   baryon number   PDG encoding
    //                          stable         lifetime     decay table
    //                      shortlived          subType   anti_encoding
    //                      excitation 
    anInstance = new G4Ions(      name,     2991.17*MeV, 2.501e-12*MeV,   +1.0*eplus,
                                     1,              +1,             0,
                                     0,               0,             0,
                             "nucleus",               0,            +3,   1010010030,
                                 false,       0.2631*ns,       nullptr,
		                 false,        "static",   -1010010030, 
		                   0.0,                0                             );
    // Magnetic Moment
    G4double mN = eplus*hbar_Planck/2.0/(proton_mass_c2 /c_squared);
    anInstance->SetPDGMagneticMoment( 2.97896248 * mN );

    // create Decay Table 
    G4DecayTable* table = new G4DecayTable;
    // create decay channels
    const G4double half_br_lambda_to_p_pim = 0.5*0.639;
    const G4double half_br_lambda_to_n_piz = 0.5*0.358;
    G4VDecayChannel** mode = new G4VDecayChannel*[4];
    // lambda -> proton + pi- , with 50% probability of capturing the proton
    mode[0] = new G4PhaseSpaceDecayChannel( "hypertriton", half_br_lambda_to_p_pim, 3,
					    "deuteron", "proton", "pi-" );
    mode[1] = new G4PhaseSpaceDecayChannel( "hypertriton", half_br_lambda_to_p_pim, 2,
					    "He3", "pi-" );
    // lambda -> neutron + pi0 , with 50% probability of capturing the neutron
    mode[2] = new G4PhaseSpaceDecayChannel( "hypertriton", half_br_lambda_to_n_piz, 3,
					    "deuteron", "neutron", "pi0" );
    mode[3] = new G4PhaseSpaceDecayChannel( "hypertriton", half_br_lambda_to_n_piz, 2,
					    "triton", "pi0" );
    for ( G4int index = 0; index < 4; ++index ) table->Insert( mode[index] );
    delete [] mode;
    anInstance->SetDecayTable( table );
  }
  theInstance = static_cast< G4HyperTriton* >( anInstance );
  return theInstance;
}


G4HyperTriton* G4HyperTriton::HyperTritonDefinition() {
  return Definition();
}


G4HyperTriton* G4HyperTriton::HyperTriton() {
  return Definition();
}
