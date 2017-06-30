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
// $Id: G4EtaPrime.cc 102905 2017-03-02 09:50:56Z gcosmo $
//
// 
// ----------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, 8 June 1998 Hisaya Kurashige
// **********************************************************************
//  New impelemenataion as an utility class  M.Asai, 26 July 2004
// ----------------------------------------------------------------------
//      Update mass (PDG2006)              Oct. 11 2006 H.Kurashige 
//

#include "G4EtaPrime.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                         EtaPrime                               ###
// ######################################################################

G4EtaPrime* G4EtaPrime::theInstance = 0;

G4EtaPrime* G4EtaPrime::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "eta_prime";
  // search in particle table]
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance ==0)
  {
  // create particle
  //
  //    Arguments for constructor are as follows
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table
  //             shortlived      subType    anti_encoding

   anInstance = new G4ParticleDefinition(
                 name,     0.95778*GeV,     0.197*MeV,         0.0,
                    0,              -1,            +1,
                    0,               0,            +1,
              "meson",               0,             0,         331,
                false,          0.0*ns,          NULL,
                false,     "eta_prime",           331);
 //create Decay Table
  G4DecayTable* table = new G4DecayTable();

 // create decay channels
  G4VDecayChannel** mode = new G4VDecayChannel*[5];
  // EtaPrime -> eta + pi+ + pi-
  mode[0] = new G4PhaseSpaceDecayChannel("eta_prime",0.429,3,"eta","pi+","pi-");
  // EtaPrime -> eta + pi0 + pi0
  mode[1] = new G4PhaseSpaceDecayChannel("eta_prime",0.222,3,"eta","pi0","pi0");
  // EtaPrime -> rho0 + gamma
  mode[2] = new G4PhaseSpaceDecayChannel("eta_prime",0.291,2,"rho0","gamma");

  // EtaPrime -> gamma + gamma
  mode[3] = new G4PhaseSpaceDecayChannel("eta_prime",0.0220,2,"gamma","gamma");

  // EtaPrime -> omega + gamma
  mode[4] = new G4PhaseSpaceDecayChannel("eta_prime",0.0275,2,"omega","gamma");

  for (G4int index=0; index <5; index++ ) table->Insert(mode[index]);
  delete [] mode;

   anInstance->SetDecayTable(table);
  }
  theInstance = reinterpret_cast<G4EtaPrime*>(anInstance);
  return theInstance;
}

G4EtaPrime*  G4EtaPrime::EtaPrimeDefinition()
{
  return Definition();
}

G4EtaPrime*  G4EtaPrime::EtaPrime()
{
  return Definition();
}

