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
// $Id: G4Globals.hh,v 1.4 2001-07-11 09:59:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
