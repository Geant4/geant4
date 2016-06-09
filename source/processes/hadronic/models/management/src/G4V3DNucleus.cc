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
// $Id: G4V3DNucleus.cc,v 1.4 2005/06/04 13:40:04 jwellisc Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
#include "G4V3DNucleus.hh"
#include "G4HadronicException.hh"

G4V3DNucleus::G4V3DNucleus()
{
}

G4V3DNucleus::G4V3DNucleus(const G4V3DNucleus &)
{
}


G4V3DNucleus::~G4V3DNucleus()
{
}


const G4V3DNucleus & G4V3DNucleus::operator=(const G4V3DNucleus &)
{
  G4String text = "G4V3DNucleus::operator= meant to not be accessable";
  throw G4HadronicException(__FILE__, __LINE__, text); 
  return *this;
}


int G4V3DNucleus::operator==(const G4V3DNucleus &) const
{
  return 0;
}

int G4V3DNucleus::operator!=(const G4V3DNucleus &) const
{
  return 1;
}


