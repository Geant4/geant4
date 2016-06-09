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
// $Id: G4SigmacPlusPlus.cc,v 1.10 2005/01/14 03:49:11 asaim Exp $
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

#include "G4SigmacPlusPlus.hh"
#include "G4ParticleTable.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                      SigmacPlusPlus                            ###
// ######################################################################

G4SigmacPlusPlus* G4SigmacPlusPlus::theInstance = 0;

G4SigmacPlusPlus* G4SigmacPlusPlus::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "sigma_c++";
  // search in particle table]
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (theInstance ==0)
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
                 name,    2.4526*GeV,       2.0*MeV,  +2.0*eplus,
                    1,              +1,             0,
                    2,              +2,             0,
             "baryon",               0,            +1,        4222,
                false,          0.0*ns,          NULL,
                false,       "sigma_c");
 //create Decay Table
  G4DecayTable* table = new G4DecayTable();

  // create decay channels
  G4VDecayChannel** mode = new G4VDecayChannel*[1];
  // sigma_c++ -> lambda_c+ + pi+
  mode[0] = new G4PhaseSpaceDecayChannel("sigma_c++",1.000,2,"lambda_c+","pi+");

  for (G4int index=0; index <1; index++ ) table->Insert(mode[index]);
  delete [] mode;

   anInstance->SetDecayTable(table);
  }
  theInstance = reinterpret_cast<G4SigmacPlusPlus*>(anInstance);
  return theInstance;
}

G4SigmacPlusPlus*  G4SigmacPlusPlus::SigmacPlusPlusDefinition()
{
  return Definition();
}

G4SigmacPlusPlus*  G4SigmacPlusPlus::SigmacPlusPlus()
{
  return Definition();
}


