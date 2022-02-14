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

#include "G4AntiHyperH4.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

// ###############################################################
// ###                   ANTI  HYPERH4                         ###
// ###############################################################

G4AntiHyperH4* G4AntiHyperH4::theInstance = nullptr;


G4AntiHyperH4* G4AntiHyperH4::Definition() {
  if ( theInstance != nullptr ) return theInstance;
  const G4String name = "anti_hyperH4";
  // search in particle table
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4Ions* anInstance =  reinterpret_cast< G4Ions* >( pTable->FindParticle( name ) );
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
    anInstance = new G4Ions(      name,     3922.40*MeV, 2.501e-12*MeV,   -1.0*eplus,
                                     4,              +1,             0,
                                     0,               0,             0,
                        "anti_nucleus",               0,            -4,  -1010010040,
                                 false,       0.2631*ns,       nullptr,
		                 false,        "static",    1010010040, 
		                   0.0,                0                             );
    // Magnetic Moment
    G4double mN = eplus*hbar_Planck/2.0/(proton_mass_c2 /c_squared);
    anInstance->SetPDGMagneticMoment( 2.97896248 * mN );

    // create Decay Table 
    G4DecayTable* table = new G4DecayTable;
    // create decay channels
    /*
    // The decay "mode[3]" produces the secondary "anti_H4" which is not existing 
    // in Geant4: we therefore skip it for the time being (similarly for hyperH4,
    // to keep the symmetry between particle and anti-particle).
    const G4double half_br_lambda_to_p_pim = 0.5*0.639;
    const G4double half_br_lambda_to_n_piz = 0.5*0.358;
    G4VDecayChannel** mode = new G4VDecayChannel*[4];
    // anti_lambda -> anti_proton + pi+ , with 50% probability of capturing the anti_proton
    mode[0] = new G4PhaseSpaceDecayChannel( "anti_hyperH4", half_br_lambda_to_p_pim, 3,
                                            "anti_triton", "anti_proton", "pi+" );
    mode[1] = new G4PhaseSpaceDecayChannel( "anti_hyperH4", half_br_lambda_to_p_pim, 2,
                                            "anti_alpha", "pi+" );
    // anti_lambda -> anti_neutron + pi0 , with 50% probability of capturing the anti_neutron
    mode[2] = new G4PhaseSpaceDecayChannel( "anti_hyperH4", half_br_lambda_to_n_piz, 3,
                                            "anti_triton", "anti_neutron", "pi0" );
    mode[3] = new G4PhaseSpaceDecayChannel( "anti_hyperH4", half_br_lambda_to_n_piz, 2,
                                            "anti_H4", "pi0" );
    for ( G4int index = 0; index < 4; ++index ) table->Insert( mode[index] );
    */
    // Replacement decay for the time being    
    const G4double half_br_lambda_to_p_pim = 0.5*0.639;
    const G4double br_lambda_to_n_piz = 0.358;
    G4VDecayChannel** mode = new G4VDecayChannel*[3];
    // anti_lambda -> anti_proton + pi+ , with 50% probability of capturing the anti_proton
    mode[0] = new G4PhaseSpaceDecayChannel( "anti_hyperH4", half_br_lambda_to_p_pim, 3,
                                            "anti_triton", "anti_proton", "pi+" );
    mode[1] = new G4PhaseSpaceDecayChannel( "anti_hyperH4", half_br_lambda_to_p_pim, 2,
                                            "anti_alpha", "pi+" );
    // anti_lambda -> anti_neutron + pi0 , with 0% probability of capturing the anti_neutron
    mode[2] = new G4PhaseSpaceDecayChannel( "anti_hyperH4", br_lambda_to_n_piz, 3,
                                            "anti_triton", "anti_neutron", "pi0" );
    for ( G4int index = 0; index < 3; ++index ) table->Insert( mode[index] );
    //---
    delete [] mode;
    anInstance->SetDecayTable( table );
  }
  theInstance = reinterpret_cast< G4AntiHyperH4* >( anInstance );
  return theInstance;
}


G4AntiHyperH4* G4AntiHyperH4::AntiHyperH4Definition() {
  return Definition();
}


G4AntiHyperH4* G4AntiHyperH4::AntiHyperH4() {
  return Definition();
}
