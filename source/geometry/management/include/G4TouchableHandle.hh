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
// G4TouchableHandle
//
// Class description:
//
// A type providing reference counting mechanism for any kind of touchable
// objects.
// The basic rule for the use of this type is that this handle must always
// be exchanged by reference and never be dinamically allocated (i.e. never
// instantiated using 'new').
// For more details see G4ReferenceCountedHandle.

// Author: Radovan Chytracek (CERN), 22.10.2001
// ----------------------------------------------------------------------
#ifndef G4TOUCHABLEHANDLE_HH
#define G4TOUCHABLEHANDLE_HH

#include "G4VTouchable.hh"
#include "G4ReferenceCountedHandle.hh"

/**
 * @brief G4TouchableHandle is a type providing reference counting mechanism
 * for any kind of touchable objects.
 * The basic rule for the use of this type is that the handle must always
 * be exchanged by reference and never be dinamically allocated (i.e. never
 * instantiated using 'new').
 */

using G4TouchableHandle = G4ReferenceCountedHandle<G4VTouchable>;

#endif
