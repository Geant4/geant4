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
// $Id:$
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// 
// G4MultiUnion
//
// Class description:
//
//   Typedef for UMultiUnion from USolids module.

// History:
// 13.11.13 G.Cosmo, CERN/PH
// --------------------------------------------------------------------
#ifndef G4MULTIUNION_HH
#define G4MULTIUNION_HH

#include "G4UMultiUnion.hh"

#if defined(G4GEOM_USE_USOLIDS)

typedef G4UMultiUnion G4MultiUnion;

#else

typedef void G4MultiUnion;

#endif
#endif
