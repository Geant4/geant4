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
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      4th April 1996, G.Cosmo
// **********************************************************************
//  New impelemenataion as an utility class             M.Asai, 26 July 2004
//  Updated mass, width and lifetime following PDG 2023 S.Okada, 5 November 2023
// ----------------------------------------------------------------------

#include "G4XicPlus.hh"

#include "G4ParticleTable.hh"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"

G4XicPlus* G4XicPlus::theInstance = nullptr;

G4XicPlus* G4XicPlus::Definition()
{
  if (theInstance != nullptr) return theInstance;
  const G4String name = "xi_c+";
  // search in particle table]
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == nullptr) {
    // create particle
    //
    //    Arguments for constructor are as follows
    //               name             mass          width         charge
    //             2*spin           parity  C-conjugation
    //          2*Isospin       2*Isospin3       G-parity
    //               type    lepton number  baryon number   PDG encoding
    //             stable         lifetime    decay table
    //             shortlived      subType    anti_encoding

    // clang-format off
   anInstance = new G4ParticleDefinition(
                 name,     2.46771*GeV,   1.45e-9*MeV,   +1.*eplus,
                    1,              +1,             0,
                    1,              +1,             0,
             "baryon",               0,            +1,        4232,
                false,    0.4531e-3*ns,       nullptr,
                false,       "xi_c");
    // clang-format on
  }
  theInstance = static_cast<G4XicPlus*>(anInstance);
  return theInstance;
}

G4XicPlus* G4XicPlus::XicPlusDefinition()
{
  return Definition();
}

G4XicPlus* G4XicPlus::XicPlus()
{
  return Definition();
}
