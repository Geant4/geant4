// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayTrajectoryPoint.hh,v 1.3 2000-03-09 17:38:33 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
    inline int operator==(const G4RayTrajectoryPoint& right) const
    { return (this==&right); };

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
    inline void SetSurfaceNormal(G4ThreeVector val) { surfaceNormal = val; }
    inline G4ThreeVector GetSurfaceNormal() const { return surfaceNormal; }
    inline void SetStepLength(G4double val) { stepLength = val; }
    inline G4double GetStepLength() const { return stepLength; }
  
};

extern G4Allocator<G4RayTrajectoryPoint> G4RayTrajectoryPointAllocator;

inline void* G4RayTrajectoryPoint::operator new(size_t)
{
   void *aTrajectoryPoint;
   aTrajectoryPoint = (void *) G4RayTrajectoryPointAllocator.MallocSingle();
   return aTrajectoryPoint;
}

inline void G4RayTrajectoryPoint::operator delete(void *aTrajectoryPoint)
{
   G4RayTrajectoryPointAllocator.FreeSingle((G4RayTrajectoryPoint *) aTrajectoryPoint);
}

#endif

