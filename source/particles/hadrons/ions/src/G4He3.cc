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
// $Id: G4He3.cc,v 1.9 2004-09-02 01:52:34 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      24th April 1998, H.Kurashige
// **********************************************************************
//  New impelemenataion as an utility class  M.Asai, 26 July 2004
// ----------------------------------------------------------------------

#include "G4He3.hh"
#include "G4ParticleTable.hh"

// ######################################################################
// ###                           He3                                  ###
// ######################################################################

G4ParticleDefinition* G4He3::theInstance = 0;

G4ParticleDefinition* G4He3::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "He3";
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
                 name,   2.80923*GeV,       0.0*MeV,  +2.0*eplus,
                    1,              +1,             0,
                    0,               0,             0,
            "nucleus",               0,            +3,           0,
                 true,            -1.0,          NULL,
             false,           "static"
              );

  theInstance->SetAtomicNumber(2);
  theInstance->SetAtomicMass(3);
  return theInstance;
}

G4ParticleDefinition*  G4He3::He3Definition()
{
  return Definition();
}

G4ParticleDefinition*  G4He3::He3()
{
  return Definition();
}


