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
// $Id: G4Point3DVector.hh,v 1.6 2001-07-11 09:59:36 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4Point3DVector
//
// Class description:
//
// A value collection of points in 3D space (G4Point3D).

// Authors: J.Sulkimo, P.Urban.
// ----------------------------------------------------------------------
#ifndef included_G4Point3DVector
#define included_G4Point3DVector

#include "g4std/vector"
#include "G4Point3D.hh"

typedef G4std::vector<G4Point3D> G4Point3DVector;

#endif
