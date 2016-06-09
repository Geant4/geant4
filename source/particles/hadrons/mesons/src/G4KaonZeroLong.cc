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
// $Id: G4KaonZeroLong.cc,v 1.15 2005/01/14 03:49:16 asaim Exp $
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

#include "G4KaonZeroLong.hh"
#include "G4ParticleTable.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4KL3DecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                      KAONZEROLONG                              ###
// ######################################################################

G4KaonZeroLong* G4KaonZeroLong::theInstance = 0;

G4KaonZeroLong* G4KaonZeroLong::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "kaon0L";
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
                 name,    0.497672*GeV,  1.273e-14*MeV,         0.0,
                    0,              -1,             0,
                    1,               0,             0,
              "meson",               0,             0,         130,
                false,         51.7*ns,          NULL,
                false,          "kaon",           130);

 //create Decay Table
  G4DecayTable* table = new G4DecayTable();

 // create decay channels
  G4VDecayChannel** mode = new G4VDecayChannel*[6];
  // kaon0L -> pi0 + pi0 + pi0
  mode[0] = new G4PhaseSpaceDecayChannel("kaon0L",0.211,3,"pi0","pi0","pi0");
  // kaon0L -> pi0 + pi+ + pi-
  mode[1] = new G4PhaseSpaceDecayChannel("kaon0L",0.126,3,"pi0","pi+","pi-");
  // kaon0L -> pi- + e+ + nu_e (Ke3)
  mode[2] = new G4KL3DecayChannel("kaon0L",0.194,"pi-","e+","nu_e");
  // kaon0L -> pi+ + e- + anti_nu_e (Ke3)
  mode[3] = new G4KL3DecayChannel("kaon0L",0.194,"pi+","e-","anti_nu_e");
  // kaon0L -> pi- + mu+ + nu_mu (Kmu3)
  mode[4] = new G4KL3DecayChannel("kaon0L",0.136,"pi-","mu+","nu_mu");
  // kaon0L -> pi+ + mu- + anti_nu_mu (Kmu3)
  mode[5] = new G4KL3DecayChannel("kaon0L",0.136,"pi+","mu-","anti_nu_mu");

  for (G4int index=0; index <6; index++ ) table->Insert(mode[index]);
  delete [] mode;

   anInstance->SetDecayTable(table);
  }
  theInstance = reinterpret_cast<G4KaonZeroLong*>(anInstance);
  return theInstance;
}

G4KaonZeroLong*  G4KaonZeroLong::KaonZeroLongDefinition()
{
  return Definition();
}

G4KaonZeroLong*  G4KaonZeroLong::KaonZeroLong()
{
  return Definition();
}

