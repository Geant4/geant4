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
// Thread Local Storage typedefs

// Author: G.Cosmo, 01.10.2012 - Created
// --------------------------------------------------------------------

// Fundamental definitions
#ifndef G4GMAKE
#  include "G4GlobalConfig.hh"
#endif

#ifndef G4_TLS
#  define G4_TLS 1

#  if defined(G4MULTITHREADED)
#    if(defined(__MACH__) && defined(__clang__)) ||                            \
      (defined(__linux__) && defined(__clang__))
#      define G4ThreadLocalStatic static thread_local
#      define G4ThreadLocal thread_local
#    elif((defined(__linux__) || defined(__MACH__)) &&                         \
          !defined(__INTEL_COMPILER) && defined(__GNUC__) &&                   \
          (__GNUC__ >= 4 && __GNUC_MINOR__ < 9))
#      define G4ThreadLocalStatic static __thread
#      define G4ThreadLocal thread_local
#    elif((defined(__linux__) || defined(__MACH__)) &&                         \
            !defined(__INTEL_COMPILER) && defined(__GNUC__) &&                 \
            (__GNUC__ >= 4 && __GNUC_MINOR__ >= 9) ||                          \
          __GNUC__ >= 5)
#      define G4ThreadLocalStatic static thread_local
#      define G4ThreadLocal thread_local
#    elif((defined(__linux__) || defined(__MACH__)) &&                         \
          defined(__INTEL_COMPILER))
#      if __INTEL_COMPILER >= 1500
#        define G4ThreadLocalStatic static thread_local
#        define G4ThreadLocal thread_local
#      else
#        define G4ThreadLocalStatic static __thread
#        define G4ThreadLocal __thread
#      endif
#    elif defined(_AIX)
#      define G4ThreadLocalStatic static thread_local
#      define G4ThreadLocal thread_local
#    elif defined(WIN32)
#      define G4ThreadLocalStatic static thread_local
#      define G4ThreadLocal thread_local
#    else
#      error                                                                   \
        "No Thread Local Storage (TLS) technology supported for this platform. Use sequential build !"
#    endif
#  else
#    define G4ThreadLocalStatic static
#    define G4ThreadLocal
#  endif

#endif
