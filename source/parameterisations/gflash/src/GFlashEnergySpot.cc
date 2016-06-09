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
// $Id: GFlashEnergySpot.cc,v 1.4 2005/10/04 09:08:33 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
// ------------------------------------------------------------
// GEANT 4 class implementation
//
//      ---------------- GFlashEnergySpot ----------------
//
// Author: Joanna Weng - 9.11.2004
// ------------------------------------------------------------

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Polyline.hh"
#include "G4VVisManager.hh"
#include "G4Step.hh"
#include "GFlashEnergySpot.hh"

GFlashEnergySpot::GFlashEnergySpot() {}

GFlashEnergySpot::GFlashEnergySpot(const G4ThreeVector& point, G4double E)
{
  Point = point;
  Energy = E;
  // initialize shower start @@@@@
}

GFlashEnergySpot::~GFlashEnergySpot() {}
