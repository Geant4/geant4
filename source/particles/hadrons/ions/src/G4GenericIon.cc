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
// $Id: G4GenericIon.cc,v 1.10 2005/01/14 03:49:13 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      4 Dec. 98 by H.Kurashige
// ----------------------------------------------------------------------
//  New impelemenataion as an utility class  M.Asai, 26 July 2004
// ----------------------------------------------------------------------

#include "G4GenericIon.hh"
#include "G4ParticleTable.hh"

// ######################################################################
// ###                           GenericIon                           ###
// ######################################################################
G4GenericIon* G4GenericIon::theInstance = 0;

G4GenericIon* G4GenericIon::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "GenericIon";
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
//!!!!
//!!!! this particle should not be used for tracking
//!!!! all properties except name and type are meaningless
//!!!!
   anInstance = new G4ParticleDefinition(
                 name,   0.9382723*GeV,       0.0*MeV,       eplus,
                    1,              +1,             0,          
                    1,              +1,             0,             
	    "nucleus",               0,            +1,           0,
		 true,            -1.0,          NULL,
             false,           "generic"
              );
  }
  theInstance = reinterpret_cast<G4GenericIon*>(anInstance);
  return theInstance;
}

G4GenericIon*  G4GenericIon::GenericIonDefinition()
{ 
  return Definition();
}

G4GenericIon*  G4GenericIon::GenericIon()
{ 
  return Definition();
}

