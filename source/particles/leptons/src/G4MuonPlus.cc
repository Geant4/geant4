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
// $Id: G4MuonPlus.cc,v 1.9 2004-09-02 01:52:40 asaim Exp $
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

#include "G4MuonPlus.hh"
#include "G4ParticleTable.hh"

#include "G4MuonDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                          MUONPLUS                              ###
// ######################################################################
G4ParticleDefinition* G4MuonPlus::theInstance = 0;

G4ParticleDefinition* G4MuonPlus::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "mu+";
  // search in particle table]
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
                 name,   0.1056584*GeV, 2.99591e-16*MeV,   +1.*eplus, 
		    1,               0,             0,          
		    0,               0,             0,             
	     "lepton",              -1,             0,         -13,
		false,      2197.03*ns,          NULL,
             false,           "mu"
              );
  //create Decay Table 
  G4DecayTable* table = new G4DecayTable();
  // create a decay channel
  G4VDecayChannel* mode = new G4MuonDecayChannel("mu+",1.00);
  table->Insert(mode);
  theInstance->SetDecayTable(table);

  return theInstance;
}

G4ParticleDefinition*  G4MuonPlus::MuonPlusDefinition()
{
  return Definition();
}

G4ParticleDefinition*  G4MuonPlus::MuonPlus()
{
  return Definition();
}

