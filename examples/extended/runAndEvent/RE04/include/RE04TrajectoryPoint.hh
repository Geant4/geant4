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
/// \file runAndEvent/RE04/include/RE04TrajectoryPoint.hh
/// \brief Definition of the RE04TrajectoryPoint class
//
//
#ifndef RE04TrajectoryPoint_h
#define RE04TrajectoryPoint_h 1

#include "G4VTrajectoryPoint.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"
class G4Material;

//
/// Trajectory point class
///
/// - new, delete and "==" operators are overwritten
///
/// - const G4ThreeVector GetPosition() const
///    gets the position of this trajectory
///
/// - const G4Material* GetMaterial() const
///    gets material that this trajectory has.
///
/// - const std::map<G4String,G4AttDef>* GetAttDefs() const
///    defines the position and the material as attiributes 
///
/// - std::vector<G4AttValue>* CreateAttValues() const
///    sets and returns the attributes
//
////////////////////////
class RE04TrajectoryPoint : public G4VTrajectoryPoint
//////////////////////// 
{

//--------
public: 
//--------

// Constructor/Destructor
   RE04TrajectoryPoint();
   RE04TrajectoryPoint(G4ThreeVector pos,const G4Material* mat);
   RE04TrajectoryPoint(const RE04TrajectoryPoint &right);
   virtual ~RE04TrajectoryPoint();

// Operators
   inline void *operator new(size_t);
   inline void operator delete(void *aTrajectoryPoint);
   inline G4bool operator==(const RE04TrajectoryPoint& right) const
   { return (this==&right); };

// Get/Set functions
   inline virtual const G4ThreeVector GetPosition() const
   { return fPosition; };
   inline const G4Material* GetMaterial() const
   { return fpMaterial; };

// Get method for HEPRep style attributes
   virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
   virtual std::vector<G4AttValue>* CreateAttValues() const;

//---------
   private:
//---------

// Member data
   G4ThreeVector fPosition;
   const G4Material* fpMaterial;

};

extern G4ThreadLocal G4Allocator<RE04TrajectoryPoint> * faTrajPointAllocator;

inline void* RE04TrajectoryPoint::operator new(size_t)
{
  if(!faTrajPointAllocator)
    faTrajPointAllocator = new G4Allocator<RE04TrajectoryPoint>;
  return (void *) faTrajPointAllocator->MallocSingle();
}

inline void RE04TrajectoryPoint::operator delete(void *aTrajectoryPoint)
{
   faTrajPointAllocator->FreeSingle((RE04TrajectoryPoint *) aTrajectoryPoint);
}

#endif

