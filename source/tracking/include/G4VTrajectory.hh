// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VTrajectory.hh,v 1.1 1999-04-15 08:47:10 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4VTrajectory.hh
//
// Description:
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
   public:

   G4VTrajectory() {;}
   virtual ~G4VTrajectory() {;}

   inline int operator == (const G4VTrajectory& right){return (this==&right);}

   virtual void ShowTrajectory() const = 0;
   virtual void DrawTrajectory(G4int i_mode=0) const = 0;
   virtual void AppendStep(const G4Step*) = 0;
   virtual int GetPointEntries() const = 0;
   virtual G4VTrajectoryPoint* GetPoint(int) const = 0;
};

#endif










