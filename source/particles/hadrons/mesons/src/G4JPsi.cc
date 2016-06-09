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
// $Id: G4JPsi.cc,v 1.14 2005/01/14 03:49:16 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      Created                 Hisaya Kurashige, 16 June 1997
// **********************************************************************
//  New impelemenataion as an utility class  M.Asai, 26 July 2004
// ----------------------------------------------------------------------

#include "G4JPsi.hh"
#include "G4ParticleTable.hh"

// ######################################################################
// ###                                JPsi                            ###
// ######################################################################

G4JPsi* G4JPsi::theInstance = 0;

G4JPsi* G4JPsi::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "J/psi";
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
                 name,    3.09688*GeV,     0.087*MeV,          0.,
                    2,              -1,            -1,
                    0,               0,            -1,
              "meson",               0,             0,         443,
                false,          0.0*ns,          NULL,
                false,         "J/psi",           443);
  }
  theInstance = reinterpret_cast<G4JPsi*>(anInstance);
  return theInstance;
}

G4JPsi*  G4JPsi::JPsiDefinition()
{
  return Definition();
}

G4JPsi*  G4JPsi::JPsi()
{
  return Definition();
}


