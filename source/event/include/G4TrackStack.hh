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
// $Id: G4TrackStack.hh,v 1.5 2001-07-13 15:01:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  Last Modification : 09/Dec/96 M.Asai
//


#ifndef G4TrackStack_h
#define G4TrackStack_h 1

#include "G4StackedTrack.hh"
#include "globals.hh"

// class description:
//
//  This is a stack class used by G4StackManager. This class object
// stores G4StackedTrack class objects in the form of bi-directional
// linked list.

class G4TrackStack 
{
  public:
      G4TrackStack();
      ~G4TrackStack();

  private:
      const G4TrackStack & operator=
                          (const G4TrackStack &right);
      G4int operator==(const G4TrackStack &right) const;
      G4int operator!=(const G4TrackStack &right) const;

  public:
      void PushToStack(G4StackedTrack * aStackedTrack);
      G4StackedTrack * PopFromStack();
      void GrabFromStack(G4StackedTrack * aStackedTrack);
      void clear();
      void TransferTo(G4TrackStack * aStack);

  private:
      G4int n_stackedTrack;
      G4StackedTrack * firstStackedTrack;
      G4StackedTrack * lastStackedTrack;

  public:
      inline G4int GetNTrack() const
      { return n_stackedTrack; }
};

#endif

