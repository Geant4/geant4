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
// $Id: G4AntiOmegaMinus.cc,v 1.11 2004-09-02 01:52:31 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#include "G4AntiOmegaMinus.hh"
#include "G4ParticleTable.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                           AntiOmegaMinus                       ###
// ######################################################################

G4ParticleDefinition* G4AntiOmegaMinus::theInstance = 0;

G4ParticleDefinition* G4AntiOmegaMinus::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "anti_omega-";
  // search in particle table
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  theInstance = pTable->FindParticle(name);
  if (theInstance !=0) return theInstance;

  // create particle
  //
  //    Arguments for constructor are as follows
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity   
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table
  //             shortlived      subType    anti_encoding
  
  theInstance = new G4ParticleDefinition(
                 name,    1.67245*GeV,  8.02e-12*MeV,       eplus,
                    3,              +1,             0,
                    0,               0,             0,
             "baryon",               0,            -1,       -3334,
                false,       0.0822*ns,          NULL,
                false,       "omega");
 //create Decay Table 
  G4DecayTable* table = new G4DecayTable();
  
  // create decay channels
  G4VDecayChannel** mode = new G4VDecayChannel*[3];
  // anti_omega- -> anti_lambda + kaon+
  mode[0] = new G4PhaseSpaceDecayChannel("anti_omega-",0.678,2,"anti_lambda","kaon+");
  // anti_omega- -> anti_xi0 + pi+
  mode[1] = new G4PhaseSpaceDecayChannel("anti_omega-",0.236,2,"anti_xi0","pi+");
  // anti_omega- -> anti_xi- + pi0
  mode[2] = new G4PhaseSpaceDecayChannel("anti_omega-",0.086,2,"anti_xi-","pi0");

  for (G4int index=0; index <3; index++ ) table->Insert(mode[index]);
  delete [] mode;

  theInstance->SetDecayTable(table);
  return theInstance;
}

G4ParticleDefinition*  G4AntiOmegaMinus::AntiOmegaMinusDefinition()
{ 
  return Definition();
}

G4ParticleDefinition*  G4AntiOmegaMinus::AntiOmegaMinus()
{ 
  return Definition();
}

