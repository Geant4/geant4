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
//      Update width (PDG2006)              Oct. 11 2006 H.Kurashige 
//

#include "G4KaonZeroLong.hh"
#include "G4SystemOfUnits.hh"
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
                 name,    0.497614*GeV,  1.287e-14*MeV,         0.0,
                    0,              -1,             0,
                    1,               0,             0,
              "meson",               0,             0,         130,
                false,        51.16*ns,          NULL,
                false,          "kaon",           130);

 //create Decay Table
  G4DecayTable* table = new G4DecayTable();

 // create decay channels
  G4VDecayChannel** mode = new G4VDecayChannel*[6];
  // kaon0L -> pi0 + pi0 + pi0
  mode[0] = new G4PhaseSpaceDecayChannel("kaon0L",0.1952,3,"pi0","pi0","pi0");
  // kaon0L -> pi0 + pi+ + pi-
  mode[1] = new G4PhaseSpaceDecayChannel("kaon0L",0.1254,3,"pi0","pi+","pi-");
  // kaon0L -> pi- + e+ + nu_e (Ke3)
  mode[2] = new G4KL3DecayChannel("kaon0L",0.2027,"pi-","e+","nu_e");
  // kaon0L -> pi+ + e- + anti_nu_e (Ke3)
  mode[3] = new G4KL3DecayChannel("kaon0L",0.2027,"pi+","e-","anti_nu_e");
  // kaon0L -> pi- + mu+ + nu_mu (Kmu3)
  mode[4] = new G4KL3DecayChannel("kaon0L",0.1352,"pi-","mu+","nu_mu");
  // kaon0L -> pi+ + mu- + anti_nu_mu (Kmu3)
  mode[5] = new G4KL3DecayChannel("kaon0L",0.1352,"pi+","mu-","anti_nu_mu");

  for (G4int index=0; index <6; index++ ) table->Insert(mode[index]);
  delete [] mode;

   anInstance->SetDecayTable(table);
  }
  theInstance = static_cast<G4KaonZeroLong*>(anInstance);
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

