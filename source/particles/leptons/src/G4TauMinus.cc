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
// $Id: G4TauMinus.cc,v 1.11 2004-09-02 01:52:40 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      4th April 1996, G.Cosmo
//      7 July 1996                   H.Kurashige 
// **********************************************************************
//  New impelemenataion as an utility class  M.Asai, 26 July 2004
// ----------------------------------------------------------------------

#include "G4TauMinus.hh"
#include "G4ParticleTable.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4TauLeptonicDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                          TAUMINUS                             ###
// ######################################################################
G4ParticleDefinition* G4TauMinus::theInstance = 0;

G4ParticleDefinition* G4TauMinus::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "tau-";
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
                 name,     1.77705*GeV,  2.265e-9*MeV,    -1.*eplus, 
		    1,               0,             0,          
		    0,               0,             0,             
	     "lepton",               1,             0,          15,
		false,     295.6e-6*ns,          NULL,
             false,           "tau"
              );

  //create Decay Table 
  G4DecayTable* table = new G4DecayTable();

  // create decay channels
  G4VDecayChannel* mode;
  // tau- -> mu- + anti_nu_mu + nu_tau
  mode = new G4TauLeptonicDecayChannel("tau-",0.174,"mu-");
  table->Insert(mode);
  // tau- -> e- + anti_nu_e + nu_tau
  mode = new G4TauLeptonicDecayChannel("tau-",0.178,"e-");
  table->Insert(mode);
  // tau- -> pi- + nu_tau
  mode = new G4PhaseSpaceDecayChannel("tau-",0.113,2,"pi-","nu_tau");
  table->Insert(mode);
  // tau- -> pi0 + pi0 + pi- + nu_tau
  mode = new G4PhaseSpaceDecayChannel("tau-",0.252,3,"pi0","pi-","nu_tau");
  table->Insert(mode);
  // tau- -> pi0 + pi0 + pi- + nu_tau
  mode = new G4PhaseSpaceDecayChannel();
  mode->SetParent("tau-");
  mode->SetBR(0.093);
  mode->SetNumberOfDaughters(4);
  mode->SetDaughter(0,"pi0");
  mode->SetDaughter(1,"pi0");
  mode->SetDaughter(2,"pi-");
  mode->SetDaughter(3,"nu_tau");
  table->Insert(mode);
  // tau- -> pi- + pi- + pi+ + nu_tau
  mode = new G4PhaseSpaceDecayChannel();
  mode->SetParent("tau-");
  mode->SetBR(0.098);
  mode->SetNumberOfDaughters(4);
  mode->SetDaughter(0,"pi-");
  mode->SetDaughter(1,"pi-");
  mode->SetDaughter(2,"pi+");
  mode->SetDaughter(3,"nu_tau");
  table->Insert(mode);

  theInstance->SetDecayTable(table);

  return theInstance;
}

G4ParticleDefinition*  G4TauMinus::TauMinusDefinition()
{
  return Definition();
}

G4ParticleDefinition*  G4TauMinus::TauMinus()
{
  return Definition();
}


