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
// $Id: G4KaonPlus.cc,v 1.11 2005/01/14 03:49:16 asaim Exp $
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

#include "G4KaonPlus.hh"
#include "G4ParticleTable.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4KL3DecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                          KAONPLUS                              ###
// ######################################################################

G4KaonPlus* G4KaonPlus::theInstance = 0;

G4KaonPlus* G4KaonPlus::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "kaon+";
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
                 name,    0.493677*GeV,  5.315e-14*MeV,    +1.*eplus,
                    0,              -1,             0,
                    1,              +1,             0,
              "meson",               0,             0,         321,
                false,       12.371*ns,          NULL,
                false,       "kaon");

 //create Decay Table
  G4DecayTable* table = new G4DecayTable();

 // create decay channels
  G4VDecayChannel** mode = new G4VDecayChannel*[6];
  // kaon+ -> mu+ + nu_mu
  mode[0] = new G4PhaseSpaceDecayChannel("kaon+",0.635,2,"mu+","nu_mu");
  // kaon+ -> pi+ + pi0
  mode[1] = new G4PhaseSpaceDecayChannel("kaon+",0.212,2,"pi+","pi0");
  // kaon+ -> pi+ + pi+ + pi-
  mode[2] = new G4PhaseSpaceDecayChannel("kaon+",0.056,3,"pi+","pi+","pi-");
  // kaon+ -> pi+ + pi0 + pi0
  mode[3] = new G4PhaseSpaceDecayChannel("kaon+",0.017,3,"pi+","pi0","pi0");
  // kaon+ -> pi0 + e+ + nu_e (Ke3)
  mode[4] = new G4KL3DecayChannel("kaon+",0.048,"pi0","e+","nu_e");
  // kaon+ -> pi0 + mu+ + nu_mu (Kmu3)
  mode[5] = new G4KL3DecayChannel("kaon+",0.032,"pi0","mu+","nu_mu");

  for (G4int index=0; index <6; index++ ) table->Insert(mode[index]);
  delete [] mode;

   anInstance->SetDecayTable(table);
  }
  theInstance = reinterpret_cast<G4KaonPlus*>(anInstance);
  return theInstance;
}

G4KaonPlus*  G4KaonPlus::KaonPlusDefinition()
{
  return Definition();
}

G4KaonPlus*  G4KaonPlus::KaonPlus()
{
  return Definition();
}

