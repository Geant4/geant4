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
// $Id: G4VHighEnergyGenerator.cc,v 1.5 2006-06-29 20:46:03 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VHighEnergyGenerator
#include "G4VHighEnergyGenerator.hh"
#include "G4HadronicException.hh"

G4VHighEnergyGenerator::G4VHighEnergyGenerator()
{
}

G4VHighEnergyGenerator::G4VHighEnergyGenerator(const G4VHighEnergyGenerator &)
{
}


G4VHighEnergyGenerator::~G4VHighEnergyGenerator()
{
}


const G4VHighEnergyGenerator & G4VHighEnergyGenerator::operator=(const G4VHighEnergyGenerator &)
{
  G4String text = "G4VHighEnergyGenerator::operator= meant to not be accessable";
  throw G4HadronicException(__FILE__, __LINE__, text); 
  return *this;
}


int G4VHighEnergyGenerator::operator==(const G4VHighEnergyGenerator &) const
{
  return 0;
}

int G4VHighEnergyGenerator::operator!=(const G4VHighEnergyGenerator &) const
{
  return 1;
}

