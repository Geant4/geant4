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
// $Id: LXeRecorderBase.hh 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/LXe/include/LXeRecorderBase.hh
/// \brief Definition of the LXeRecorderBase class
//
// LXeRecorderBase.hh
// 1-Sep-1999 Bill Seligman

// This is an abstract base class to be used with Geant 4.0.1 (and
// possibly higher, if the User classes don't change).

// The concept of a Recorder object is that it records the activities of
// Geant in a manner that is useful to a physicist.  Perhaps this record
// takes the form of histograms, or ntuples, or entries in an Objectivity
// database.  This class does not care HOW the information is recorded; it
// abstracts the behavior of a generalized recorder of Geant variables.

// No object should be instantiated from the Recorder class (in fact, any such
// object won't do anything).  The user must define a new class (say, a class
// that creates histograms) and overload the methods of this class.

// Why do this?  First of all, it keeps all record-keeping in a single class:
// the class that inherits Recorder.  The original Geant documentation suggests
// that recording activities should be split among many different classes
// (initialization in G4UserRunAction, recording in G4UserSteppingAction, etc.).
// If you use a Recorder class, than all the record-keeping details are kept in
// a single class instead of being spread out among many different classes.

// Secondly, by using an abstract Recorder class, you hide the implementation
// details from the rest of Geant.  If you change a couple of histograms, only
// the Recorder-derived class and main() re-compile.  No other class knows or
// cares what or how you record.

// The only time this class (i.e., this header file) changes is if new
// user action classes are added to Geant.

#ifndef RECORDER_BASE_H_
#define RECORDER_BASE_H_

// The following objects are the arguments to the methods
// invoked in the user action classes.  In other words, they
// contain the variables that we are normally able to record
// in Geant.

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4Step.hh"

class LXeRecorderBase {

  public:

    virtual ~LXeRecorderBase() {};

    // The following a list of methods that correspond to the available
    // user action classes in Geant 4.0.1.   In this base class, the
    // methods are defined to do nothing.

    virtual void RecordBeginOfRun(const G4Run*) = 0;
    virtual void RecordEndOfRun(const G4Run*) = 0;
    virtual void RecordBeginOfEvent(const G4Event*) {};
    virtual void RecordEndOfEvent(const G4Event*) {};
    virtual void RecordTrack(const G4Track*) {};
    virtual void RecordStep(const G4Step*) {};

};

#endif
