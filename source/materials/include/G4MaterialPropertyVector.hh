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
//
// 
////////////////////////////////////////////////////////////////////////
//
// class G4MaterialPropertyVector
//
// Class description:
//
// A one-to-one mapping from Photon Energy to some optical property 

// File:        G4MaterialPropertyVector.hh
//
// Version:     2.0
// Created:     1996-02-08
// Author:      Juliet Armstrong
// Updated:     2011-10-13 by Peter Gumplinger
//              remove the class: simply typedef to G4PhysicsFreeVector
//
////////////////////////////////////////////////////////////////////////

#ifndef G4MaterialPropertyVector_h   
#define G4MaterialPropertyVector_h 1

/////////////
// Includes
/////////////

#include "G4PhysicsFreeVector.hh"

/////////////////////
// Class Definition
/////////////////////

using G4MaterialPropertyVector = G4PhysicsFreeVector;

#endif /* G4MaterialPropertyVector_h */
