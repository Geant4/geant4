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
// $Id: test1.cc,v 1.3 2001-07-11 10:01:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  19th July 1996

// Simply tests linkability of kernel without visualization.

#include "G4ios.hh"
#include <math.h>

#include "BuildCalorimeter.hh"
#include "G4VPhysicalVolume.hh"
#include "G4GeometryManager.hh"

#include "g4templates.hh"

main () {

  // Build some geometry and get relevant pointers.
  G4VPhysicalVolume *pWorld = BuildCalorimeter ();
  G4GeometryManager::GetInstance () -> CloseGeometry ();

  return 0;
}
