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

#include "G4KaonPlus.hh"
#include "G4SystemOfUnits.hh"
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
                 name,    0.493677*GeV, 5.317e-14*MeV,    +1.*eplus,
                    0,              -1,             0,
                    1,              +1,             0,
              "meson",               0,             0,         321,
                false,       12.380*ns,          NULL,
                false,       "kaon");

 //create Decay Table
  G4DecayTable* table = new G4DecayTable();

 // create decay channels
  G4VDecayChannel** mode = new G4VDecayChannel*[6];
  // kaon+ -> mu+ + nu_mu
  mode[0] = new G4PhaseSpaceDecayChannel("kaon+",0.6355,2,"mu+","nu_mu");
  // kaon+ -> pi+ + pi0
  mode[1] = new G4PhaseSpaceDecayChannel("kaon+",0.2066,2,"pi+","pi0");
  // kaon+ -> pi+ + pi+ + pi-
  mode[2] = new G4PhaseSpaceDecayChannel("kaon+",0.0559,3,"pi+","pi+","pi-");
  // kaon+ -> pi+ + pi0 + pi0
  mode[3] = new G4PhaseSpaceDecayChannel("kaon+",0.01761,3,"pi+","pi0","pi0");
  // kaon+ -> pi0 + e+ + nu_e (Ke3)
  mode[4] = new G4KL3DecayChannel("kaon+",0.0507,"pi0","e+","nu_e");
  // kaon+ -> pi0 + mu+ + nu_mu (Kmu3)
  mode[5] = new G4KL3DecayChannel("kaon+",0.0335,"pi0","mu+","nu_mu");

  for (G4int index=0; index <6; index++ ) table->Insert(mode[index]);
  delete [] mode;

   anInstance->SetDecayTable(table);
  }
  theInstance = static_cast<G4KaonPlus*>(anInstance);
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

