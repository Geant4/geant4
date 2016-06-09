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
// $Id: G4Version.hh,v 1.1 2005/11/28 10:59:30 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Version information
//
// History:
// 26.09.05 K.Murakami  - Created
//

#ifndef G4VERSION_HH
#define G4VERSION_HH

// Numbering rule for "G4VERSION_NUMBER":
// - The number is consecutive (i.e. 711) as an integer.
// - The meaning of each digit is as follows;
//
//   711
//   |--> major version number 
//    |--> minor version number
//     |--> patch number

#ifndef G4VERSION_NUMBER
#define G4VERSION_NUMBER  800
#endif

#ifndef G4VERSION_TAG
#define G4VERSION_TAG "$Name: geant4-08-00 $"
#endif

// as variables

#include "G4String.hh"

static const G4String G4Version = "$Name: geant4-08-00 $";
static const G4String G4Date    = "(16-December-2005)";

#endif
