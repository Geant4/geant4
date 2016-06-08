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
// $Id: G4PVHit.cc,v 1.4.8.1 2001/06/28 19:11:25 gunter Exp $
// GEANT4 tag $Name:  $
//

#include "G4PVHit.hh"
#include "globals.hh"

G4PVHit::G4PVHit()
{;}

G4PVHit::~G4PVHit()
{;}

int G4PVHit::operator==(const G4PVHit &right) const
{ return false; }

void G4PVHit::Draw()
{;}

void G4PVHit::Print()
{;}

