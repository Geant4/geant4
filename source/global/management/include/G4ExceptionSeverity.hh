//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ExceptionSeverity.hh,v 1.2 2002-08-26 14:56:46 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class Description:
//
// Specifies the severity of G4Exception
//
//  FatalException
//   Error is severe or it happens at the initialization time.
//   Program should be aborted and core dump will be generated.
//
//  FatalErrorInArgument
//   Fatal error caused by most likely the mis-use of interfaces
//   by the user's code. Program should be aborted and core dump 
//   will be generated.
//
//  RunMustBeAborted
//   Error happens at initialization of a run (ex. at the 
//   moment of closing geometry), or some unpleasant situation
//   occurs during the event loop. Current run will be aborted
//   and the application returns to "Idle" state.
//
//  EventMustBeAborted
//   Error happens during tracking a particle. The event currently
//   being processed should be aborted, run will not be aborted.
//
//  JustWarning
//   Just display messages.
//

#ifndef G4ExceptionSeverity_H
#define G4ExceptionSeverity_H 1

enum G4ExceptionSeverity 
  { FatalException,
    FatalErrorInArgument,
    RunMustBeAborted,
    EventMustBeAborted,
    JustWarning };

#endif

