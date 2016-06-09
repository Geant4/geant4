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
// $Id: G4VHighEnergyGenerator.cc,v 1.3 2003/11/03 17:54:18 hpw Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
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

