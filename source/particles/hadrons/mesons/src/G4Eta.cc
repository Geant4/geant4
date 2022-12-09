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
//      Update mass (PDG2006)              Oct. 11 2006 H.Kurashige 
//

#include "G4Eta.hh"
#include "G4SystemOfUnits.hh"
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
                 name,    0.547862*GeV,      1.31*keV,         0.0,
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
  mode[0] = new G4PhaseSpaceDecayChannel("eta",0.3942,2,"gamma","gamma");
  // eta -> pi0 + pi0 + pi0
  mode[1] = new G4PhaseSpaceDecayChannel("eta",0.3256,3,"pi0","pi0","pi0");
  // eta -> pi0 + pi+ + pi-
  mode[2] = new G4PhaseSpaceDecayChannel("eta",0.226,3,"pi0","pi+","pi-");
  // eta -> gamma + pi+ + pi-
  mode[3] = new G4PhaseSpaceDecayChannel("eta",0.0468,3,"gamma","pi+","pi-");

  for (G4int index=0; index <4; index++ ) table->Insert(mode[index]);
  delete [] mode;

   anInstance->SetDecayTable(table);
  }
  theInstance = static_cast<G4Eta*>(anInstance);
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

