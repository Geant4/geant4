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
// $Id: $
// GEANT4 tag $Name: $
//
//
// Thread Local Storage typedefs

// History:
// 01.10.2012 G.Cosmo - Created
 
#ifndef G4_TLS
#define G4_TLS

#if defined (G4LIB_BUILD_MULTITHREAD)
  #if defined(__MACH__) && defined(__clang__) && defined(MAC_OS_X_VERSION_10_7) && defined(__x86_64__)
  #  define G4TLOCAL static __thread
  #elif defined(__linux__) || defined(_AIX)
  #  define G4TLOCAL static __thread
  #elif defined(WIN32)
  #  define G4TLOCAL static __declspec(thread)
  #else
  #  error "No Thread Local Storage (TLS) technology supported for this platform. Use sequential build !"
  #endif
#else
  #  define G4TLOCAL static
#endif

#endif
