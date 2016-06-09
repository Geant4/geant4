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
// $Id: G4V3DNucleus.cc,v 1.5 2002/12/12 19:17:30 gunter Exp $
// GEANT4 tag $Name: geant4-05-01 $
//
#include "G4V3DNucleus.hh"

G4V3DNucleus::G4V3DNucleus()
{
}

G4V3DNucleus::G4V3DNucleus(const G4V3DNucleus &right)
{
}


G4V3DNucleus::~G4V3DNucleus()
{
}


const G4V3DNucleus & G4V3DNucleus::operator=(const G4V3DNucleus &right)
{
  G4Exception("G4V3DNucleus::operator= meant to not be accessable"); // needs to be looked at @@
  return *this;
}


int G4V3DNucleus::operator==(const G4V3DNucleus &right) const
{
  return 0;
}

int G4V3DNucleus::operator!=(const G4V3DNucleus &right) const
{
  return 1;
}


