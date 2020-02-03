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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4AtomicShellEnumerator.hh
//
// Author:        Alfonso Mantero 
//
// Creation date: 15.03.2011
//
// Modifications:
//
// Class Description:
//
// Definition table between atomic shell names and EADL Shell Ids
// -------------------------------------------------------------------
//

#ifndef G4AtomicShellEnumerator_h
#define G4AtomicShellEnumerator_h 1

enum G4AtomicShellEnumerator // not following EADL syntax
{
  //***********EADL shell Id              
  fKShell = 0,  //  1
  fL1Shell = 1, //  3
  fL2Shell = 2,	//  5
  fL3Shell = 3,	//  6
  fM1Shell = 4,	//  8
  fM2Shell = 5,	//  10
  fM3Shell = 6,	//  11
  fM4Shell = 7,	//  13
  fM5Shell = 8  //  14
};

#endif
