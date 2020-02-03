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
// ****************************************************************
//  New impelemenataion as an utility class  H.Kurashige, 14 July 2004
// ----------------------------------------------------------------

#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"


// ######################################################################
// ###                          GEANTINO                              ###
// ######################################################################
G4Geantino* G4Geantino::theInstance = 0;

G4Geantino*  G4Geantino::Definition() 
{
  if (theInstance !=0) return theInstance;

  const G4String name = "geantino";
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
		  name,           0.0*MeV,       0.0*MeV,         0.0, 
		  0,                    0,             0,          
		  0,                    0,             0,             
		  "geantino",           0,             0,           0,
		  true,              -1.0,          NULL,
		  false,        "geantino",            0
		);
  }
  theInstance = reinterpret_cast<G4Geantino*>(anInstance);
  return theInstance;

}

G4Geantino*  G4Geantino::GeantinoDefinition() 
{
  return Definition();
}

G4Geantino*  G4Geantino::Geantino() 
{
  return Definition();
}
