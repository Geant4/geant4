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
// $Id: G4Deuteron.cc,v 1.10 2005/01/14 03:49:13 asaim Exp $
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

#include "G4Deuteron.hh"
#include "G4ParticleTable.hh"

// ######################################################################
// ###                           DEUTERON                             ###
// ######################################################################

G4Deuteron* G4Deuteron::theInstance = 0;

G4Deuteron* G4Deuteron::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "deuteron";
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
                 name,   1.875613*GeV,       0.0*MeV,  +1.0*eplus,
                    2,              +1,             0,
                    0,               0,             0,
            "nucleus",               0,            +2,           0,
                 true,            -1.0,          NULL,
             false,           "static"
              );

   anInstance->SetAtomicNumber(1);
   anInstance->SetAtomicMass(2);
  }
  theInstance = reinterpret_cast<G4Deuteron*>(anInstance);
  return theInstance;
}

G4Deuteron*  G4Deuteron::DeuteronDefinition()
{
  return Definition();
}

G4Deuteron*  G4Deuteron::Deuteron()
{
  return Definition();
}


