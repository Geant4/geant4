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
// ----------------------------------------------------------------------
// G4TiMemory
//
// Provides empty macros when Geant4 is compiled with TiMemory disabled
// ----------------------------------------------------------------------

#ifndef g4timemory_hh_
#define g4timemory_hh_

#include "globals.hh"

//----------------------------------------------------------------------------//
#ifdef GEANT4_USE_TIMEMORY

#   if defined __GNUC__
#       pragma GCC diagnostic push
#       pragma GCC diagnostic ignored "-Wexceptions"
#       pragma GCC diagnostic ignored "-Wunused-private-field"
#   endif

#include <timemory/timemory.hpp>

typedef tim::auto_timer G4AutoTimer;

inline void InitializeTiMemory()
{
    tim::manager* instance = tim::manager::instance();
    instance->enable(true);
}

#   if defined __GNUC__
#       pragma GCC diagnostic pop
#   endif

#else

#define TIMEMORY_AUTO_TIMER(str)
#define TIMEMORY_AUTO_TIMER_OBJ(str) {}

#define TIMEMORY_BASIC_AUTO_TIMER(str)
#define TIMEMORY_BASIC_AUTO_TIMER_OBJ(str) {}

#define TIMEMORY_DEBUG_BASIC_AUTO_TIMER(str)
#define TIMEMORY_DEBUG_AUTO_TIMER(str)

inline void InitializeTiMemory()
{ }

#endif
//----------------------------------------------------------------------------//

#endif

