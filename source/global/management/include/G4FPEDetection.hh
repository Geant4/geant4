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
// $Id: G4FPEDetection.hh,v 1.3 2010-07-14 10:45:46 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include <stdlib.h>  // abort(), exit()

static void TerminationSignalHandler(int sig)
{
  std::cerr << "ERROR: " << sig;
  std::string message;
  switch (SIGFPE)
  {
  #ifdef FPE_NOOP  // occurs in OS X
    case FPE_NOOP:
      message = "No Operation.";
      break;
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
  std::cerr << " - " << message << std::endl;

  ::abort();
}

#ifdef __linux__
#ifdef __GNUC__
  #include <features.h>
  #include <fenv.h>
  #include <csignal>

  struct sigaction termaction, oldaction;

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

    sigset_t *def_set;
    def_set=&termaction.sa_mask;
    sigfillset(def_set);
    sigdelset(def_set,SIGFPE);
    termaction.sa_handler=TerminationSignalHandler;
    termaction.sa_flags=0;
    sigaction(SIGFPE, &termaction,&oldaction);
  }
#endif
#elif __MACH__      // MacOSX

  #include <fenv.h>
  #include <signal.h>

  #define DEFINED_PPC      (defined(__ppc__) || defined(__ppc64__))
  #define DEFINED_INTEL    (defined(__i386__) || defined(__x86_64__))

  #if DEFINED_PPC

    #define FE_EXCEPT_SHIFT 22  // shift flags right to get masks
    #define FM_ALL_EXCEPT    FE_ALL_EXCEPT >> FE_EXCEPT_SHIFT 

    static int feenableexcept (unsigned int excepts)
    {
      static fenv_t fenv;
      unsigned int new_excepts = (excepts & FE_ALL_EXCEPT) >> FE_EXCEPT_SHIFT,
                   old_excepts;  // all previous masks

      if ( fegetenv (&fenv) )  { return -1; }
      old_excepts = (fenv & FM_ALL_EXCEPT) << FE_EXCEPT_SHIFT;
      fenv = (fenv & ~new_excepts) | new_excepts;

      return ( fesetenv (&fenv) ? -1 : old_excepts );
    }

    static int fedisableexcept (unsigned int excepts)
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

    static int feenableexcept (unsigned int excepts)
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

    static int fedisableexcept (unsigned int excepts)
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

  #endif  // PPC or INTEL enabling

  static void InvalidOperationDetection()
  {
    struct sigaction termaction;

    std::cout << std::endl
              << "        "
              << "############################################" << std::endl
              << "        "
              << "!!! WARNING - FPE detection is activated !!!" << std::endl
              << "        "
              << "############################################" << std::endl;

    feenableexcept ( FE_DIVBYZERO );
    feenableexcept ( FE_INVALID   );
    fedisableexcept( FE_OVERFLOW  );
    fedisableexcept( FE_UNDERFLOW );

//    sigfillset (&termaction.sa_mask);
//    sigemptyset(&termaction.sa_mask);
    sigdelset  (&termaction.sa_mask,SIGFPE);
    termaction.sa_handler = TerminationSignalHandler;
//    termaction.sa_flags=SA_SIGINFO;
    
    sigaction(SIGFPE, &termaction, 0);
  }
#else  // Not Linux, nor MacOSX ...

  static void InvalidOperationDetection() {;}

#endif

#endif
