//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4EtaPrime.cc,v 1.16 2005/01/14 03:49:16 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// ----------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, 8 June 1998 Hisaya Kurashige
// **********************************************************************
//  New impelemenataion as an utility class  M.Asai, 26 July 2004
// ----------------------------------------------------------------------

#include "G4EtaPrime.hh"
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
                 name,    0.95777*GeV,     0.202*MeV,         0.0,
                    0,              -1,            +1,
                    0,               0,            +1,
              "meson",               0,             0,         331,
                false,          0.0*ns,          NULL,
                false,     "eta_prime",           331);
 //create Decay Table
  G4DecayTable* table = new G4DecayTable();

 // create decay channels
  G4VDecayChannel** mode = new G4VDecayChannel*[3];
  // EtaPrime -> eta + pi+ + pi-
  mode[0] = new G4PhaseSpaceDecayChannel("eta_prime",0.437,3,"eta","pi+","pi-");
  // EtaPrime -> eta + pi0 + pi0
  mode[1] = new G4PhaseSpaceDecayChannel("eta_prime",0.208,3,"eta","pi0","pi0");
  // EtaPrime -> rho0 + gamma
  mode[2] = new G4PhaseSpaceDecayChannel("eta_prime",0.302,2,"rho0","gamma");

  for (G4int index=0; index <3; index++ ) table->Insert(mode[index]);
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

