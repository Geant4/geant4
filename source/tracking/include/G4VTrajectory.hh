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
// $Id: G4VTrajectory.hh,v 1.6 2001-07-11 10:08:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4VTrajectory.hh
//
// class description:
//   This is the abstract base class of a trajectory.
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Makoto  Asai   (e-mail: asai@kekvax.kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
// ---------------------------------------------------------------

#ifndef G4VTrajectory_h
#define G4VTrajectory_h 1

class G4Step;
#include "G4VTrajectoryPoint.hh"
#include "globals.hh"

class G4VTrajectory
{
public: // without description

   G4VTrajectory() {;}
   virtual ~G4VTrajectory() {;}

   inline int operator == (const G4VTrajectory& right){return (this==&right);}

   virtual void ShowTrajectory() const = 0;
   virtual void DrawTrajectory(G4int i_mode=0) const = 0;
   virtual void AppendStep(const G4Step*) = 0;
   virtual int GetPointEntries() const = 0;
   virtual G4VTrajectoryPoint* GetPoint(int) const = 0;
   virtual void MergeTrajectory(G4VTrajectory*) = 0;
};

#endif










