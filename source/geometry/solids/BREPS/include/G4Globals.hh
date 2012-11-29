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
// $Id$
//
// Author: A.Breakstone
// Adaptation: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
//
// Defines global constants required by some surfaces
//
// ----------------------------------------------------------------------

#ifndef __G4BREPGLOBALS_H
#define __G4BREPGLOBALS_H

//
//  Define a C preprocessor constant for a Scale factor to apply to 
//  various dimensionless tests in the geometry routines which test
//  if a point is on a surface.  This number is an effective thickness
//  of a surface divided by a relevant dimension, such as the radius of
//  a cylinder.  The default value is 0.0001.
#define SURFACE_PRECISION 0.0001 
//
//  Define a C preprocessor constant for the maximum number of turns
//  allowed for a Helix, which is used in some of the geometry routines
//  to limit the size of for or while loops.  The default value is 50.
#define HELIX_MAX_TURNS 50

#endif
