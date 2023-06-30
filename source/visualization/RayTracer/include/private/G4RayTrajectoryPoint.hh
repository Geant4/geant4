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
//
//

// class description:
//
//  This is a concrete class of G4VTrajectoryPoint which represents a point where
// a ray crosses upon a surface of a volume regardless of its visibility. Objects
// of this class are created by G4RayTrajectory class object.
//

/////////////////////////
//G4RayTrajectoryPoint.hh
/////////////////////////

#ifndef G4RayTrajectoryPoint_h
#define G4RayTrajectoryPoint_h 1

class G4VisAttributes;
#include "globals.hh"
#include "G4VTrajectoryPoint.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class G4RayTrajectoryPoint :public G4VTrajectoryPoint
{
  public:
    G4RayTrajectoryPoint();
    virtual ~G4RayTrajectoryPoint();

    inline void *operator new(size_t);
    inline void operator delete(void *aTrajectoryPoint);
  //    inline G4bool operator==(const G4RayTrajectoryPoint& right) const
  // { return (this==&right); };

  private:
    const G4VisAttributes* preStepAtt;
    const G4VisAttributes* postStepAtt;
    G4ThreeVector    surfaceNormal;
    G4double         stepLength;

  public:
    inline void SetPreStepAtt(const G4VisAttributes* val) { preStepAtt = val; }
    inline const G4VisAttributes* GetPreStepAtt() const { return preStepAtt; }
    inline void SetPostStepAtt(const G4VisAttributes* val) { postStepAtt = val; }
    inline const G4VisAttributes* GetPostStepAtt() const { return postStepAtt; }
    inline void SetSurfaceNormal(const G4ThreeVector& val) { surfaceNormal = val; }
    inline G4ThreeVector GetSurfaceNormal() const { return surfaceNormal; }
    inline void SetStepLength(G4double val) { stepLength = val; }
    inline G4double GetStepLength() const { return stepLength; }

    inline const G4ThreeVector GetPosition() const { return G4ThreeVector();}
    // Dummy function (not used) to satisfy base class pure virtual function.
};

#if defined G4RAYTRACER_ALLOC_EXPORT
  extern G4DLLEXPORT G4Allocator<G4RayTrajectoryPoint>*& rayTrajectoryPointAllocator();
#else
  extern G4DLLIMPORT G4Allocator<G4RayTrajectoryPoint>*& rayTrajectoryPointAllocator();
#endif

inline void* G4RayTrajectoryPoint::operator new(size_t)
{
   if(!rayTrajectoryPointAllocator())
   { rayTrajectoryPointAllocator() = new G4Allocator<G4RayTrajectoryPoint>; }
   return (void *) rayTrajectoryPointAllocator()->MallocSingle();
}

inline void G4RayTrajectoryPoint::operator delete(void *aTrajectoryPoint)
{
   rayTrajectoryPointAllocator()->FreeSingle((G4RayTrajectoryPoint *) aTrajectoryPoint);
}

#endif

