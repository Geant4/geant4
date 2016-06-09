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
// INCL++ revision: v5.0_rc3
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#ifndef G4INCLFPEDebug_hh
#define G4INCLFPEDebug_hh 1

#ifdef INCL_FPE_DEBUG

#ifdef __linux__
#ifdef __GNUC__
#include <features.h>
#include <fenv.h>
#include <csignal>
#include <cstdlib>
#include <iostream>

struct sigaction termaction, oldaction;

void TerminationSignalHandler(G4int sig)
{
  std::cerr << "ERROR: " << sig;
  std::string message;
  switch (SIGFPE)
    {
    case FPE_INTDIV:
      message = "Integer divide by zero.";
      break;
    case FPE_INTOVF:
      message = "Integer overflow.";
      break;
    case FPE_FLTDIV:
      message = "Floating poG4int divide by zero.";
      break;
    case FPE_FLTOVF:
      message = "Floating poG4int overflow.";
      break;
    case FPE_FLTUND:
      message = "Floating poG4int underflow.";
      break;
    case FPE_FLTRES:
      message = "Floating poG4int inexact result.";
      break;
    case FPE_FLTINV:
      message = "Floating poG4int invalid operation.";
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

void enableFPEDetection ()
{
  std::cerr << std::endl
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
#else
void enableFPEDetection() {
  std::cerr <<"FPE detection is not supported by your platform." << std::endl;
}
#endif

#endif

#endif
