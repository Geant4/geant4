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
// $Id: meshdefs.hh,v 1.4 2001-07-11 09:59:19 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Tube/Cone Meshing constants for extent calculations

// History:
// 13.08.95 P.Kent Created separate file

#ifndef MESHDEFS_HH
#define MESHDEFS_HH

#include "globals.hh"

const G4double kMeshAngleDefault=(M_PI/4); // Angle for mesh `wedges' in rads
                                 // Works best when simple fraction of M_PI/2

const G4int kMinMeshSections=3;	 // Min wedges+1 to make
const G4int kMaxMeshSections=37; // max wedges+1 to make
                                 // =>10 degrees/wedge for complete tube

#endif
