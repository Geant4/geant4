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
// $Id: G4RayTrajectory.hh,v 1.16 2010-05-29 21:09:40 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

// class description:
//
//  This is a concrete class of G4VTrajectory which represents a trajectory of a ray.
//  This class is used by G4RayTracer. Objects of this class are created by G4RTTrackingAction.
//  This class does not have concrete implementations of draw and print methods but
// information stored in this class is used by determining a colour of a pixel of the picture.
//

///////////////////
//G4RayTrajectory.hh
///////////////////

#ifndef G4RayTrajectory_h
#define G4RayTrajectory_h 1

class G4Step;

#include "G4VTrajectory.hh"
#include "G4Allocator.hh"
#include <stdlib.h>
#include <vector>
#include "globals.hh"
#include "G4Track.hh"
#include "G4RayTrajectoryPoint.hh"


class G4RayTrajectory : public G4VTrajectory
{
   public:

   G4RayTrajectory(); 
   G4RayTrajectory(G4RayTrajectory & right);
   virtual ~G4RayTrajectory();

   inline void* operator new(size_t);
   inline void  operator delete(void*);
  //   inline int operator == (const G4RayTrajectory& right){return (this==&right);}

   virtual void AppendStep(const G4Step*);
   virtual void ShowTrajectory(std::ostream&) const;
   virtual void DrawTrajectory() const {;}
   virtual void DrawTrajectory(G4int) const {;}
   virtual int GetPointEntries() const {return positionRecord->size();}
   virtual G4VTrajectoryPoint* GetPoint(G4int i) const 
   { return (*positionRecord)[i]; }
   G4RayTrajectoryPoint* GetPointC(G4int i) const 
   { return (*positionRecord)[i]; }
   virtual void MergeTrajectory(G4VTrajectory* secondTrajectory);  

// Get/Set functions to satisfy pure virtual functions of base class.
   inline G4int GetTrackID() const { return 0; }
   inline G4int GetParentID() const { return 0; }
   inline G4String GetParticleName() const { return ""; }
   inline G4double GetCharge() const { return 0.; }
   inline G4int GetPDGEncoding() const { return 0; }
   inline G4ThreeVector GetInitialMomentum() const { return G4ThreeVector(); }

   private:

   std::vector<G4RayTrajectoryPoint*>* positionRecord;
};

#if defined G4VIS_ALLOC_EXPORT
  extern G4DLLEXPORT G4Allocator<G4RayTrajectory> G4RayTrajectoryAllocator;
#else
  extern G4DLLIMPORT G4Allocator<G4RayTrajectory> G4RayTrajectoryAllocator;
#endif

inline void* G4RayTrajectory::operator new(size_t)
{
   void* aTrajectory;
   aTrajectory = (void*)G4RayTrajectoryAllocator.MallocSingle();
   return aTrajectory;
}

inline void G4RayTrajectory::operator delete(void* aTrajectory)
{
   G4RayTrajectoryAllocator.FreeSingle((G4RayTrajectory*)aTrajectory);
}


#endif

