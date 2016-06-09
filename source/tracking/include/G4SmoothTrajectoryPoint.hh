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
// $Id: G4SmoothTrajectoryPoint.hh,v 1.15 2010-10-27 07:57:21 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

