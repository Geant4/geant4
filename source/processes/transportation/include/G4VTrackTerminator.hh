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
// $Id: G4VTrackTerminator.hh,v 1.1 2002-08-29 15:32:37 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4VTrackTerminator
//
// Class description:
//
// This is an interface for an object which maybe told to kill a track.
// The type it provides is needed in case importance biasing and
// scoring is done at the same time. 
// For navigation in the parallel geometry which is done 
// by the importance biasing process (it derives from G4ParallelTransport)
// importance biasing has to be done before scoring (to score the correct
// volume). But since scoring would not be called if the importance
// biasing kills a track it only tells a G4VTrackTerminator to kill the 
// track. The scoring process implements the interface G4VTrackTerminator
// and is now able to kill the track after it has done the scoring.
// 
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4VTrackTerminator_hh
#define G4VTrackTerminator_hh G4VTrackTerminator_hh

class G4VTrackTerminator {
public:
  virtual ~G4VTrackTerminator(){}
  virtual void KillTrack() = 0;
};

#endif
