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
// $Id: G4ExceptionSeverity.hh,v 1.1 2002-08-19 18:20:11 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4ExceptionSeverity_H
#define G4ExceptionSeverity_H 1

// Class Description:
//
// Specifies the severity of G4Exception
//
//  FatalException
//   Error is so sevire or it happens at the initialization time.
//   Program should be aborted and core dump is to be generated.
//
//  FatalErrorInArgument
//   Fatal error caused by most likely the miss-use of the toolkit
//   by the user's code. Program should be aborted and core dump 
//   is to be generated.
//
//  RunMustBeAborted
//   Error happens at the initialization of a run, i.e. at the 
//   moment of closing geometry. Or some unpleasant situation
//   occurs during the event loop. Current run will be aborted
//   and program goes back to Idle state.
//
//  EventMustBeAborted
//   Error happens during tracking a particle. Currently processing
//   event should be aborted, but remaining event could be processed.
//
//  JustWarning
//   Just display messages.
//

enum G4ExceptionSeverity 
  {FatalException, FatalErrorInArgument, RunMustBeAborted,
                                 EventMustBeAborted, JustWarning};

#endif

