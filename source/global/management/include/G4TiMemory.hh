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
// G4TiMemory
//
// Description:
//
// Provides empty macros when Geant4 is compiled with TiMemory disabled

// Author: Jonathan R. Madsen, 25 April 2019
// --------------------------------------------------------------------
#ifndef g4timemory_hh
#define g4timemory_hh 1

// Fundamental definitions
#ifndef G4GMAKE
#  include "G4GlobalConfig.hh"
#endif

#include "globals.hh"

//----------------------------------------------------------------------------//
#ifdef GEANT4_USE_TIMEMORY

#  include <timemory/timemory.hpp>

using G4AutoTimer = tim::auto_timer;

#else

#  include <ostream>
#  include <string>

namespace tim
{
  template <typename... _Args>
  void timemory_init(_Args...)
  {}
  inline void timemory_finalize() {}
  inline void print_env() {}

  /// this provides "functionality" for *_HANDLE macros
  /// and can be omitted if these macros are not utilized
  struct dummy
  {
    template <typename... _Args>
    dummy(_Args&&...)
    {}
    ~dummy()            = default;
    dummy(const dummy&) = default;
    dummy(dummy&&)      = default;
    dummy& operator=(const dummy&) = default;
    dummy& operator=(dummy&&) = default;

    void start() {}
    void stop() {}
    void conditional_start() {}
    void conditional_stop() {}
    void report_at_exit(bool) {}
    template <typename... _Args>
    void mark_begin(_Args&&...)
    {}
    template <typename... _Args>
    void mark_end(_Args&&...)
    {}
    friend std::ostream& operator<<(std::ostream& os, const dummy&)
    {
      return os;
    }
  };

}  // namespace tim

// startup/shutdown/configure
#  define TIMEMORY_INIT(...)
#  define TIMEMORY_FINALIZE()
#  define TIMEMORY_CONFIGURE(...)

// label creation
#  define TIMEMORY_BASIC_LABEL(...) std::string("")
#  define TIMEMORY_LABEL(...) std::string("")
#  define TIMEMORY_JOIN(...) std::string("")

// define an object
#  define TIMEMORY_BLANK_MARKER(...)
#  define TIMEMORY_BASIC_MARKER(...)
#  define TIMEMORY_MARKER(...)

// define an unique pointer object
#  define TIMEMORY_BLANK_POINTER(...)
#  define TIMEMORY_BASIC_POINTER(...)
#  define TIMEMORY_POINTER(...)

// define an object with a caliper reference
#  define TIMEMORY_BLANK_CALIPER(...)
#  define TIMEMORY_BASIC_CALIPER(...)
#  define TIMEMORY_CALIPER(...)

// define a static object with a caliper reference
#  define TIMEMORY_STATIC_BLANK_CALIPER(...)
#  define TIMEMORY_STATIC_BASIC_CALIPER(...)
#  define TIMEMORY_STATIC_CALIPER(...)

// invoke member function on caliper reference or type within reference
#  define TIMEMORY_CALIPER_APPLY(...)
#  define TIMEMORY_CALIPER_TYPE_APPLY(...)

// get an object
#  define TIMEMORY_BLANK_HANDLE(...) tim::dummy()
#  define TIMEMORY_BASIC_HANDLE(...) tim::dummy()
#  define TIMEMORY_HANDLE(...) tim::dummy()

// get a pointer to an object
#  define TIMEMORY_BLANK_POINTER_HANDLE(...) nullptr
#  define TIMEMORY_BASIC_POINTER_HANDLE(...) nullptr
#  define TIMEMORY_POINTER_HANDLE(...) nullptr

// debug only
#  define TIMEMORY_DEBUG_BLANK_MARKER(...)
#  define TIMEMORY_DEBUG_BASIC_MARKER(...)
#  define TIMEMORY_DEBUG_MARKER(...)

// auto-timers
#  define TIMEMORY_BLANK_AUTO_TIMER(...)
#  define TIMEMORY_BASIC_AUTO_TIMER(...)
#  define TIMEMORY_AUTO_TIMER(...)
#  define TIMEMORY_BLANK_AUTO_TIMER_HANDLE(...)
#  define TIMEMORY_BASIC_AUTO_TIMER_HANDLE(...)
#  define TIMEMORY_AUTO_TIMER_HANDLE(...)
#  define TIMEMORY_DEBUG_BASIC_AUTO_TIMER(...)
#  define TIMEMORY_DEBUG_AUTO_TIMER(...)

using G4AutoTimer = tim::dummy;

#endif

//----------------------------------------------------------------------------//

#endif
