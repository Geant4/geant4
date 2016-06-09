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
// $Id: G4ios.hh,v 1.9 2004/06/09 07:30:01 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 
// ---------------------------------------------------------------
// GEANT 4 class header file
//
// G4ios.hh
//
// ---------------------------------------------------------------
#ifndef included_G4ios
#define included_G4ios

#include "G4Types.hh"

#include <iostream>

#if defined G4IOS_EXPORT
  extern G4DLLEXPORT std::ostream G4cout;
  extern G4DLLEXPORT std::ostream G4cerr;
#else
  extern G4DLLIMPORT std::ostream G4cout;
  extern G4DLLIMPORT std::ostream G4cerr;
#endif

#define G4cin std::cin
#define G4endl std::endl

#endif
