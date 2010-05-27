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
// $Id: ParametrisedBox.cc,v 1.5 2010-05-27 15:00:18 allison Exp $
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
 const G4VPhysicalVolume*) const
{
  box.SetXHalfLength ( 0.5 * m * (n + 1));
  box.SetYHalfLength ( 1 * m * (n + 1));
  box.SetZHalfLength ( 1.5 * m * (n + 1));
}
