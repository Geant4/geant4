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
// $Id: ParametrisedBox.cc,v 1.3 2001-07-11 10:09:28 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Test of parametrisation.  Based on Hans-Peter Wellisch's test.

#include "ParametrisedBox.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"

void ParametrisedBox::ComputeTransformation
(const G4int n,
 G4VPhysicalVolume* pRep) const
{
  pRep->SetTranslation (G4ThreeVector (0, n * 2 * m, 0));
}

void ParametrisedBox::ComputeDimensions
(G4Box& box,
 const G4int n,
 const G4VPhysicalVolume* pRep) const
{
  box.SetXHalfLength ( 0.5 * m * (n + 1));
  box.SetYHalfLength ( 1 * m * (n + 1));
  box.SetZHalfLength ( 1.5 * m * (n + 1));
}
