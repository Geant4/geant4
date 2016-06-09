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
// $Id: G4PionZero.cc,v 1.15 2005/01/14 03:49:16 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      4th April 1996, G.Cosmo
// **********************************************************************
//  New impelemenataion as an utility class  M.Asai, 26 July 2004
// ----------------------------------------------------------------------

#include "G4PionZero.hh"
#include "G4ParticleTable.hh"


#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DalitzDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                          PIONZERO                              ###
// ######################################################################

G4PionZero* G4PionZero::theInstance = 0;

G4PionZero* G4PionZero::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "pi0";
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
                 name,    0.1349764*GeV,   7.8e-06*MeV,         0.0,
                    0,              -1,            +1,
                    2,               0,            -1,
              "meson",               0,             0,         111,
                false,       8.4e-8*ns,          NULL,
                false,            "pi",          111);

 //create Decay Table
  G4DecayTable* table = new G4DecayTable();

  // create a decay channel
  G4VDecayChannel* mode;
  // pi0 -> gamma + gamma
  mode = new G4PhaseSpaceDecayChannel("pi0",0.988,2,"gamma","gamma");
  table->Insert(mode);
  // pi0 -> gamma + e+ + e-
  mode = new G4DalitzDecayChannel("pi0",0.012,"e-","e+");
  table->Insert(mode);

   anInstance->SetDecayTable(table);
  }
  theInstance = reinterpret_cast<G4PionZero*>(anInstance);
  return theInstance;
}

G4PionZero*  G4PionZero::PionZeroDefinition()
{
  return Definition();
}

G4PionZero*  G4PionZero::PionZero()
{
  return Definition();
}

