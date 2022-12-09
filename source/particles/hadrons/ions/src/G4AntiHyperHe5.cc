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

#include "G4AntiHyperHe5.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

// ##################################################################
// ###                   ANTI  HYPERHE5                           ###
// ##################################################################

G4AntiHyperHe5* G4AntiHyperHe5::theInstance = nullptr;


G4AntiHyperHe5* G4AntiHyperHe5::Definition() {
  if ( theInstance != nullptr ) return theInstance;
  const G4String name = "anti_hyperHe5";
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
    anInstance = new G4Ions(      name,     4840.00*MeV, 2.501e-12*MeV,   -2.0*eplus,
                                     3,              +1,             0,
                                     0,               0,             0,
                        "anti_nucleus",               0,            -5,  -1010020050,
                                 false,       0.2631*ns,       nullptr,
		                 false,        "static",    1010020050, 
		                   0.0,                0                             );
    // Magnetic Moment
    G4double mN = eplus*hbar_Planck/2.0/(proton_mass_c2 /c_squared);
    anInstance->SetPDGMagneticMoment( 2.97896248 * mN );

    // create Decay Table
    G4DecayTable* table = new G4DecayTable;
    // create decay channels
    /*
    // The decays "mode[1]" and "mode[3]" produce the secondaries "anti_Li5" and "anti_He5",
    // respectively, which are not existing in Geant4: we therefore skip them for the time being
    // (similarly for hyperHe5, to keep the symmetry between particle and anti-particle).
    const G4double half_br_lambda_to_p_pim = 0.5*0.639;
    const G4double half_br_lambda_to_n_piz = 0.5*0.358;
    G4VDecayChannel** mode = new G4VDecayChannel*[4];
    // anti_lambda -> anti_proton + pi+ , with 50% probability of capturing the anti_proton
    mode[0] = new G4PhaseSpaceDecayChannel( "anti_hyperHe5", half_br_lambda_to_p_pim, 3,
                                            "anti_alpha", "anti_proton", "pi+" );
    mode[1] = new G4PhaseSpaceDecayChannel( "anti_hyperHe5", half_br_lambda_to_p_pim, 2,
                                            "anti_Li5", "pi+" );
    // anti_lambda -> anti_neutron + pi0 , with 50% probability of capturing the anti_neutron
    mode[2] = new G4PhaseSpaceDecayChannel( "anti_hyperHe5", half_br_lambda_to_n_piz, 3,
                                            "anti_alpha", "anti_neutron", "pi0" );
    mode[3] = new G4PhaseSpaceDecayChannel( "anti_hyperHe5", half_br_lambda_to_n_piz, 2,
                                            "anti_He5", "pi0" );
    for ( G4int index = 0; index < 4; ++index ) table->Insert( mode[index] );
    */
    // Replacement decay for the time being    
    const G4double br_lambda_to_p_pim = 0.639;
    const G4double br_lambda_to_n_piz = 0.358;
    G4VDecayChannel** mode = new G4VDecayChannel*[2];
    // anti_lambda -> anti_proton + pi+ , with 0% probability of capturing the anti_proton
    mode[0] = new G4PhaseSpaceDecayChannel( "anti_hyperHe5", br_lambda_to_p_pim, 3,
                                            "anti_alpha", "anti_proton", "pi+" );
    // anti_lambda -> anti_neutron + pi0 , with 0% probability of capturing the anti_neutron
    mode[1] = new G4PhaseSpaceDecayChannel( "anti_hyperHe5", br_lambda_to_n_piz, 3,
                                            "anti_alpha", "anti_neutron", "pi0" );
    for ( G4int index = 0; index < 2; ++index ) table->Insert( mode[index] );    
    //---
    delete [] mode;
    anInstance->SetDecayTable( table );
  }
  theInstance = static_cast< G4AntiHyperHe5* >( anInstance );
  return theInstance;
}


G4AntiHyperHe5* G4AntiHyperHe5::AntiHyperHe5Definition() {
  return Definition();
}


G4AntiHyperHe5* G4AntiHyperHe5::AntiHyperHe5() {
  return Definition();
}
