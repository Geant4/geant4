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
// $Id: G4FPEDetection.hh 86793 2014-11-18 10:01:46Z gcosmo $
//
// 
// -*- C++ -*-
//
// -----------------------------------------------------------------------
// This global method should be used on LINUX or MacOSX platforms with gcc
// compiler for activating NaN detection and FPE signals, and forcing
// abortion of the application at the time these are detected.
// Meant to be used for debug purposes, can be activated by compiling the
// "run" module with the flag G4FPE_DEBUG set in the environment.
// -----------------------------------------------------------------------

#ifndef G4FPEDetection_h
#define G4FPEDetection_h 1

#include <iostream>
#include <stdlib.h>  /* abort(), exit() */

#if (defined(__GNUC__) && !defined(__clang__))
#ifdef __linux__
  #include <features.h>
  #include <fenv.h>
  #include <csignal>
// for G4StackBacktrace()
  #include <execinfo.h>
  #include <cxxabi.h>

  struct sigaction termaction, oldaction;

  static void G4StackBackTrace()
  {

     //   from http://linux.die.net/man/3/backtrace_symbols_fd
     #define BSIZE 50
     void * buffer[ BSIZE ];
     int nptrs = backtrace( buffer, BSIZE );
       //std::cerr << "nptrs=" << nptrs << std::endl;
     char ** strings = backtrace_symbols( buffer, nptrs );
     if ( strings == NULL ) {
       perror( "backtrace_symbols" );
       return;
     }
     std::cerr << std::endl<< "Call Stack:" << std::endl;
     for ( int j = 0; j < nptrs; j++ ){
         //printf("%s\n", strings[j]);
         std::cerr  << nptrs-j-1 <<": ";
	    //std::cerr  << strings[j] << std::endl;
       char * mangled_start = strchr( strings[j],  '(' ) + 1;
       if (mangled_start) *(mangled_start-1) = '\0'; 
       char * mangled_end   = strchr( mangled_start,'+' );
       if ( mangled_end ) *mangled_end = '\0';
         //std::cerr << "mangled .. " << mangled_start << ", len:" << strlen(mangled_start)<< std::endl;
       int status = 0;
       char *realname=0;       
       if ( mangled_end && strlen(mangled_start)  ) 
              realname = abi::__cxa_demangle( mangled_start, 0, 0, &status );
       if ( realname ) {
	 std::cerr << strings[j]<< " : " << realname  << std::endl;
	 free( realname );
       } else {
	 std::cerr << strings[j] << std::endl;
       }
     }
     free( strings );
     // c++filt can demangle: http://gcc.gnu.org/onlinedocs/libstdc++/manual/ext_demangling.html
     //-------------------------------------------------------------

  }
  static void TerminationSignalHandler(int sig, siginfo_t* sinfo, void* /* context */)
  {
    std::cerr << "ERROR: " << sig;
    std::string message = "Floating-point exception (FPE).";

    if (sinfo) {
      switch (sinfo->si_code) {
#ifdef FPE_NOOP		/* BUG: MacOSX uses this instead of INTDIV */
      case FPE_NOOP:
#endif
      case FPE_INTDIV:
        message = "Integer divide by zero.";
        break;
      case FPE_INTOVF:
        message = "Integer overflow.";
        break;
      case FPE_FLTDIV:
        message = "Floating point divide by zero.";
        break;
      case FPE_FLTOVF:
        message = "Floating point overflow.";
        break;
      case FPE_FLTUND:
        message = "Floating point underflow.";
        break;
      case FPE_FLTRES:
        message = "Floating point inexact result.";
        break;
      case FPE_FLTINV:
        message = "Floating point invalid operation.";
        break;
      case FPE_FLTSUB:
        message = "Subscript out of range.";
        break;
      default:
        message = "Unknown error.";
        break;
      }
    }

    std::cerr << " - " << message << std::endl;
    G4StackBackTrace();
    ::abort();
  }

  static void InvalidOperationDetection()
  {
    std::cout << std::endl
              << "        "
              << "############################################" << std::endl
              << "        "
              << "!!! WARNING - FPE detection is activated !!!" << std::endl
              << "        "
              << "############################################" << std::endl;

    (void) feenableexcept( FE_DIVBYZERO );
    (void) feenableexcept( FE_INVALID );
    //(void) feenableexcept( FE_OVERFLOW );
    //(void) feenableexcept( FE_UNDERFLOW );

    sigfillset(&termaction.sa_mask);
    sigdelset(&termaction.sa_mask,SIGFPE);
    termaction.sa_sigaction=TerminationSignalHandler;
    termaction.sa_flags=SA_SIGINFO;
    sigaction(SIGFPE, &termaction, &oldaction);
  }

#elif defined(__MACH__)      /* MacOSX */
  #include <fenv.h>
  #include <signal.h>

  #define DEFINED_PPC      (defined(__ppc__) || defined(__ppc64__))
  #define DEFINED_INTEL    (defined(__i386__) || defined(__x86_64__))

  #if DEFINED_PPC

    #define FE_EXCEPT_SHIFT 22  // shift flags right to get masks
    #define FM_ALL_EXCEPT    FE_ALL_EXCEPT >> FE_EXCEPT_SHIFT 

    static inline int feenableexcept (unsigned int excepts)
    {
      static fenv_t fenv;
      unsigned int new_excepts = (excepts & FE_ALL_EXCEPT) >> FE_EXCEPT_SHIFT,
                   old_excepts;  // all previous masks

      if ( fegetenv (&fenv) )  { return -1; }
      old_excepts = (fenv & FM_ALL_EXCEPT) << FE_EXCEPT_SHIFT;
      fenv = (fenv & ~new_excepts) | new_excepts;

      return ( fesetenv (&fenv) ? -1 : old_excepts );
    }

    static inline int fedisableexcept (unsigned int excepts)
    {
      static fenv_t fenv;
      unsigned int still_on = ~((excepts & FE_ALL_EXCEPT) >> FE_EXCEPT_SHIFT),
                   old_excepts;  // previous masks

      if ( fegetenv (&fenv) )  { return -1; }
      old_excepts = (fenv & FM_ALL_EXCEPT) << FE_EXCEPT_SHIFT;
      fenv &= still_on;

      return ( fesetenv (&fenv) ? -1 : old_excepts );
    }

  #elif DEFINED_INTEL

    static inline int feenableexcept (unsigned int excepts)
    {
      static fenv_t fenv;
      unsigned int new_excepts = excepts & FE_ALL_EXCEPT,
                   old_excepts;  // previous masks

      if ( fegetenv (&fenv) )  { return -1; }
      old_excepts = fenv.__control & FE_ALL_EXCEPT;

      // unmask
      //
      fenv.__control &= ~new_excepts;
      fenv.__mxcsr   &= ~(new_excepts << 7);

      return ( fesetenv (&fenv) ? -1 : old_excepts );
    }

    static inline int fedisableexcept (unsigned int excepts)
    {
      static fenv_t fenv;
      unsigned int new_excepts = excepts & FE_ALL_EXCEPT,
                   old_excepts;  // all previous masks

      if ( fegetenv (&fenv) )  { return -1; }
      old_excepts = fenv.__control & FE_ALL_EXCEPT;

      // mask
      //
      fenv.__control |= new_excepts;
      fenv.__mxcsr   |= new_excepts << 7;

      return ( fesetenv (&fenv) ? -1 : old_excepts );
    }

  #endif  /* PPC or INTEL enabling */

  static void TerminationSignalHandler(int sig, siginfo_t* sinfo, void* /* context */)
  {
    std::cerr << "ERROR: " << sig;
    std::string message = "Floating-point exception (FPE).";

    if (sinfo) {
      switch (sinfo->si_code) {
#ifdef FPE_NOOP		/* BUG: MacOSX uses this instead of INTDIV */
      case FPE_NOOP:
#endif
      case FPE_INTDIV:
        message = "Integer divide by zero.";
        break;
      case FPE_INTOVF:
        message = "Integer overflow.";
        break;
      case FPE_FLTDIV:
        message = "Floating point divide by zero.";
        break;
      case FPE_FLTOVF:
        message = "Floating point overflow.";
        break;
      case FPE_FLTUND:
        message = "Floating point underflow.";
        break;
      case FPE_FLTRES:
        message = "Floating point inexact result.";
        break;
      case FPE_FLTINV:
        message = "Floating point invalid operation.";
        break;
      case FPE_FLTSUB:
        message = "Subscript out of range.";
        break;
      default:
        message = "Unknown error.";
        break;
      }
    }

    std::cerr << " - " << message << std::endl;
    
    ::abort();
  }

  static void InvalidOperationDetection()
  {
    struct sigaction termaction, oldaction;

    std::cout << std::endl
              << "        "
              << "############################################" << std::endl
              << "        "
              << "!!! WARNING - FPE detection is activated !!!" << std::endl
              << "        "
              << "############################################" << std::endl;

    feenableexcept ( FE_DIVBYZERO );
    feenableexcept ( FE_INVALID   );
    // fedisableexcept( FE_OVERFLOW  );
    // fedisableexcept( FE_UNDERFLOW );

    sigfillset(&termaction.sa_mask);
    sigdelset(&termaction.sa_mask,SIGFPE);
    termaction.sa_sigaction=TerminationSignalHandler;
    termaction.sa_flags=SA_SIGINFO;
    sigaction(SIGFPE, &termaction, &oldaction);
  }
#else  /* Not Linux, nor MacOSX ... */

  static void InvalidOperationDetection() {;}

#endif	/* Linus or MacOSX */
#else  /* Not GCC */

  static void InvalidOperationDetection() {;}
#endif
#endif	/* G4FPEDetection_h */
