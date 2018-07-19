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
// $Id: G4EventGenerator.cc 100828 2016-11-02 15:25:59Z gcosmo $
//
// G4EventGenerator

#include "G4EventGenerator.hh"
#include "G4SystemOfUnits.hh"

G4EventGenerator::G4EventGenerator()
{
   SetMinEnergy (0 *GeV);
   SetMaxEnergy (0 *GeV);
}

G4EventGenerator::G4EventGenerator(const G4EventGenerator &) : G4HadronicInteraction()
{
}


G4EventGenerator::~G4EventGenerator()
{
}


const G4EventGenerator & G4EventGenerator::operator=(const G4EventGenerator &)
{
  throw G4HadronicException(__FILE__, __LINE__, 
                            "G4EventGenerator::operator= meant to not be accessable");
  return *this;
}


int G4EventGenerator::operator==(const G4EventGenerator &) const
{
  return 0;
}

int G4EventGenerator::operator!=(const G4EventGenerator &) const
{
  return 1;
}

