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
// INCL++ intra-nuclear cascade model
// Pekka Kaitaniemi, CEA and Helsinki Institute of Physics
// Davide Mancusi, CEA
// Alain Boudard, CEA
// Sylvie Leray, CEA
// Joseph Cugnon, University of Liege
//
// INCL++ revision: v5.1_rc11
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#ifndef G4INCLSignalHandling_hh
#define G4INCLSignalHandling_hh

#ifdef INCL_SIGNAL_HANDLING

#ifdef __linux__
#ifdef __GNUC__
#include <features.h>
#include <fenv.h>
#include <csignal>
#include <cstdlib>
#include <iostream>

G4bool INTsignalled = false, HUPsignalled = false;

struct sigaction fpeaction, hupaction, intaction, oldaction;

void FPESignalHandler(G4int /*sig*/, siginfo_t *siginfo, void * /*context*/) {
  const char *message;
  switch(siginfo->si_code)
    {
    case FPE_INTDIV:
      message = "ERROR - Integer divide by zero.          \n";
      break;
    case FPE_INTOVF:
      message = "ERROR - Integer overflow.                \n";
      break;
    case FPE_FLTDIV:
      message = "ERROR - Floating point divide by zero.   \n";
      break;
    case FPE_FLTOVF:
      message = "ERROR - Floating point overflow.         \n";
      break;
    case FPE_FLTUND:
      message = "ERROR - Floating point underflow.        \n";
      break;
    case FPE_FLTRES:
      message = "ERROR - Floating point inexact result.   \n";
      break;
    case FPE_FLTINV:
      message = "ERROR - Floating point invalid operation.\n";
      break;
    case FPE_FLTSUB:
      message = "ERROR - Subscript out of range.          \n";
      break;
    default:
      message = "ERROR - Unknown floating-point exception.\n";
      break;
    }
  ::write(STDERR_FILENO, message, 42); // the answer to the ultimate question on Life, the Universe and Everything
  ::_exit(2);
}

void INTSignalHandler(G4int /*sig*/) {
  if(INTsignalled) {
    const char *message = "Caught SIGINT while waiting for an event to finish, exiting.\n";
    ::write(STDERR_FILENO, message, 61); // the answer to the ultimate question on Life, the Universe and Everything
    ::_exit(1);
  } else {
    INTsignalled = true;
    const char *message = "Caught SIGINT, will try to quit cleanly at the end of the current event.\n";
    ::write(STDERR_FILENO, message, 74);
  }
}

void HUPSignalHandler(G4int /*sig*/) {
  HUPsignalled = true;
  const char *message = "Caught SIGHUP, will flush log and output files at the end of the current event.\n";
  ::write(STDERR_FILENO, message, 81);
}

void enableSignalHandling ()
{
  (void) feenableexcept( FE_DIVBYZERO );
  (void) feenableexcept( FE_INVALID );
  //(void) feenableexcept( FE_OVERFLOW );
  //(void) feenableexcept( FE_UNDERFLOW );

  // Set the handler for SIGFPE
  sigset_t *def_set;
  def_set=&fpeaction.sa_mask;
  sigfillset(def_set);
  sigdelset(def_set,SIGFPE);
  fpeaction.sa_sigaction=FPESignalHandler;
  fpeaction.sa_flags=SA_SIGINFO;
  sigaction(SIGFPE, &fpeaction, &oldaction);

  // Set the handler for SIGHUP
  sigset_t *def_set_hup;
  def_set_hup=&hupaction.sa_mask;
  sigfillset(def_set_hup);
  sigdelset(def_set_hup,SIGHUP);
  hupaction.sa_handler=HUPSignalHandler;
  hupaction.sa_flags=0;
  sigaction(SIGHUP, &hupaction, &oldaction);

  // Set the handler for SIGINT
  sigset_t *def_set_int;
  def_set_int=&intaction.sa_mask;
  sigfillset(def_set_int);
  sigdelset(def_set_int,SIGINT);
  intaction.sa_handler=INTSignalHandler;
  intaction.sa_flags=0;
  sigaction(SIGINT, &intaction, &oldaction);
}
#endif
#else
void enableSignalHandling() {
  std::cerr <<"Signal handling is not supported by your platform." << std::endl;
}
#endif

#endif

#endif
