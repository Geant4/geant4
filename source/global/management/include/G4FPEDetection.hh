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
// $Id: G4FPEDetection.hh,v 1.1 2006-11-15 12:09:56 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// -*- C++ -*-
//
// -----------------------------------------------------------------------
// This global method should be used on LINUX platforms with gcc compiler
// for activating NaN detection and FPE signals, and forcing abortion of
// the application at the time these are detected.
// Meant to be used for debug purposes, can be activated by compiling the
// "run" module with the flag G4FPE_DEBUG set in the environment.
// -----------------------------------------------------------------------

#ifndef G4FPEDetection_h
#define G4FPEDetection_h 1

#ifdef __GNUC__
  #include <features.h>
  #include <fenv.h>
  #include <csignal>

  #include <iostream>

  struct sigaction termaction, oldaction;

  void TerminationSignalHandler(int sig)
  {
    std::cerr << "ERROR: " << sig << std::endl;
    ::abort();
  }

  void InvalidOperationDetection()
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
#else
  void InvalidOperationDetection() {;}
#endif

#endif
