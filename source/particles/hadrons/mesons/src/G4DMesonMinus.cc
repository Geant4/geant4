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
// $Id: G4DMesonMinus.cc 102905 2017-03-02 09:50:56Z gcosmo $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//
//      Created                 Hisaya Kurashige, 16 June 1997
// **********************************************************************
//  New impelemenataion as an utility class  M.Asai, 26 July 2004
// ----------------------------------------------------------------------

#include "G4DMesonMinus.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

// ######################################################################
// ###                         DMesonMinus                            ###
// ######################################################################

G4DMesonMinus* G4DMesonMinus::theInstance = 0;

G4DMesonMinus* G4DMesonMinus::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "D-";
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
                 name,     1.86958*GeV,  6.33e-10*MeV,    -1.*eplus,
                    0,              -1,             0,
                    1,              -1,             0,
              "meson",               0,             0,        -411,
                false,     1.040e-3*ns,          NULL,
                false,       "D");
  }
  theInstance = reinterpret_cast<G4DMesonMinus*>(anInstance);
  return theInstance;
}

G4DMesonMinus*  G4DMesonMinus::DMesonMinusDefinition()
{
  return Definition();
}

G4DMesonMinus*  G4DMesonMinus::DMesonMinus()
{
  return Definition();
}

