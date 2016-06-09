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
// $Id: G4Eta.cc,v 1.15 2005/01/14 03:49:15 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      4th April 1996, G.Cosmo
//                              H.Kurashige   7 Jul 96
// **********************************************************************
//  New impelemenataion as an utility class  M.Asai, 26 July 2004
// ----------------------------------------------------------------------

#include "G4Eta.hh"
#include "G4ParticleTable.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                         ETA                                    ###
// ######################################################################

G4Eta* G4Eta::theInstance = 0;

G4Eta* G4Eta::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "eta";
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
                 name,    0.54730*GeV,      1.18*keV,         0.0,
                    0,              -1,            +1,
                    0,               0,            +1,
              "meson",               0,             0,         221,
                false,          0.0*ns,          NULL,
                false,           "eta",           221);
 //create Decay Table
  G4DecayTable* table = new G4DecayTable();

 // create decay channels
  G4VDecayChannel** mode = new G4VDecayChannel*[4];
  // eta -> gamma + gamma
  mode[0] = new G4PhaseSpaceDecayChannel("eta",0.393,2,"gamma","gamma");
  // eta -> pi0 + pi0 + pi0
  mode[1] = new G4PhaseSpaceDecayChannel("eta",0.321,3,"pi0","pi0","pi0");
  // eta -> pi0 + pi+ + pi-
  mode[2] = new G4PhaseSpaceDecayChannel("eta",0.232,3,"pi0","pi+","pi-");
  // eta -> gamma + pi+ + pi-
  mode[3] = new G4PhaseSpaceDecayChannel("eta",0.048,3,"gamma","pi+","pi-");

  for (G4int index=0; index <4; index++ ) table->Insert(mode[index]);
  delete [] mode;

   anInstance->SetDecayTable(table);
  }
  theInstance = reinterpret_cast<G4Eta*>(anInstance);
  return theInstance;
}

G4Eta*  G4Eta::EtaDefinition()
{
  return Definition();
}

G4Eta*  G4Eta::Eta()
{
  return Definition();
}

