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
// $Id: G4SmoothTrajectoryPoint.hh,v 1.10 2004/06/11 14:30:17 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
//
//---------------------------------------------------------------
//
// G4SmoothTrajectoryPoint.hh
//
// class description:
//  This class contains position and auxiliary points of a trajectory point.
//
// ---------------------------------------------------------------

#ifndef G4SmoothTrajectoryPoint_h
#define G4SmoothTrajectoryPoint_h 1

#include "G4VTrajectoryPoint.hh"
#include "globals.hh"                // Include from 'global'
#include "G4ThreeVector.hh"          // Include from 'geometry'
#include "G4Allocator.hh"            // Include from 'particle+matter'


class G4SmoothTrajectoryPoint : public G4VTrajectoryPoint
{

//--------
public: // without description
//--------

// Constructor/Destructor
   G4SmoothTrajectoryPoint();
   // No auxiliary points setter, so must set the points in the
   // constructor already (jacek 31/10/2002)
   G4SmoothTrajectoryPoint(G4ThreeVector pos,
			   std::vector<G4ThreeVector>* auxiliaryPoints);
   G4SmoothTrajectoryPoint(G4ThreeVector pos);
   G4SmoothTrajectoryPoint(const G4SmoothTrajectoryPoint &right);
   virtual ~G4SmoothTrajectoryPoint();

// Operators
   inline void *operator new(size_t);
   inline void operator delete(void *aTrajectoryPoint);
   inline int operator==(const G4SmoothTrajectoryPoint& right) const
   { return (this==&right); };

// Get/Set functions
   inline const G4ThreeVector GetPosition() const
   { return fPosition; }
   inline const std::vector<G4ThreeVector>* GetAuxiliaryPoints() const
   { return fAuxiliaryPointVector; }

// Get method for HEPRep style attributes
   virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
   virtual std::vector<G4AttValue>* CreateAttValues() const;

//---------
   private:
//---------

// Member data
   G4ThreeVector fPosition;
   std::vector<G4ThreeVector>* fAuxiliaryPointVector;
};

#if defined G4TRACKING_ALLOC_EXPORT
  extern G4DLLEXPORT G4Allocator<G4SmoothTrajectoryPoint> aSmoothTrajectoryPointAllocator;
#else
  extern G4DLLIMPORT G4Allocator<G4SmoothTrajectoryPoint> aSmoothTrajectoryPointAllocator;
#endif

inline void* G4SmoothTrajectoryPoint::operator new(size_t)
{
   void *aTrajectoryPoint;
   aTrajectoryPoint = (void *) aSmoothTrajectoryPointAllocator.MallocSingle();
   return aTrajectoryPoint;
}

inline void G4SmoothTrajectoryPoint::operator delete(void *aTrajectoryPoint)
{
   aSmoothTrajectoryPointAllocator.FreeSingle((G4SmoothTrajectoryPoint *) aTrajectoryPoint);
}

#endif

