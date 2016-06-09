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
// $Id: G4EventGenerator.cc,v 1.2 2003/11/03 17:54:53 hpw Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// G4EventGenerator
#include "G4EventGenerator.hh"

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
  throw G4HadronicException(__FILE__, __LINE__, "G4EventGenerator::operator= meant to not be accessable");
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
