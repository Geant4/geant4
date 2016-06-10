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
// ---------------------------------------------------------------
// GEANT 4 class header file
//
// Class Description:
//
// This file includes Windows declarations from <windows.h> protecting
// those defines that may cause troubles within the Geant4 code.

// ---------------------------------------------------------------
#ifndef windefs_hh
#define windefs_hh

#if defined(WIN32)
    //
	#define WIN32_LEAN_AND_MEAN
    #define NOMINMAX      // avoid redefinition of min() and max()
    #include <windows.h>
    #undef pascal         // trick to overcome redefinition of 'pascal'
    #undef scr1
    #undef scr2
    #undef rad1
    #undef rad2
    #undef small
    #undef ABSOLUTE
    #undef RELATIVE
	#undef GetObject
#endif // WIN32

#endif //windefs_hh
