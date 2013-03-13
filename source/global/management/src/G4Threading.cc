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
// $Id$
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
// G4Threading.cc
//
// ---------------------------------------------------------------
// Author: Andrea Dotti (15 Feb 2013): First Implementation
// ---------------------------------------------------------------

#include "G4Threading.hh"

#if defined(G4MULTITHREADED)
  #if defined(WIN32)
    DWORD /*WINAPI*/ G4WaitForSingleObjectInf( __in G4Mutex m ) { return WaitForSingleObject( m , INFINITE); }
    BOOL G4ReleaseMutex( __in G4Mutex m) { return ReleaseMutex(m); }
  #endif
#else
    G4int fake_mutex_lock_unlock( G4Mutex* ) { return 0; }
#endif
