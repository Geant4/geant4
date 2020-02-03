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
/// \file field/field04/include/F04TrajectoryPoint.hh
/// \brief Definition of the F04TrajectoryPoint class
//

#ifndef F04TrajectoryPoint_h
#define F04TrajectoryPoint_h 1

#include "globals.hh"

#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4TrajectoryPoint.hh"

#include "G4StepStatus.hh"

class G4Track;
class G4Step;
class G4VProcess;

class F04TrajectoryPoint : public G4TrajectoryPoint {

//--------
  public: // without description
//--------

// Constructor/Destructor

    F04TrajectoryPoint();
    F04TrajectoryPoint(const G4Track* aTrack);
    F04TrajectoryPoint(const G4Step* aStep);
    F04TrajectoryPoint(const F04TrajectoryPoint &right);
    virtual ~F04TrajectoryPoint();

// Operators

    inline void *operator new(size_t);
    inline void operator delete(void *aTrajectoryPoint);
    inline G4bool operator==(const F04TrajectoryPoint& right) const
    { return (this==&right); };

// Get/Set functions

    inline G4double GetTime() const { return fTime; };
    inline const G4ThreeVector GetMomentum() const { return fMomentum; };
    inline G4StepStatus GetStepStatus() const { return fStepStatus; };
    inline G4String GetVolumeName() const { return fVolumeName; };

// Get method for HEPRep style attributes

   virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
   virtual std::vector<G4AttValue>* CreateAttValues() const;

//---------
  private:
//---------

// Member data

    G4double      fTime;
    G4ThreeVector fMomentum;
    G4StepStatus  fStepStatus;
    G4String      fVolumeName;

};

extern G4ThreadLocal G4Allocator<F04TrajectoryPoint>* F04TrajPointAllocator;

inline void* F04TrajectoryPoint::operator new(size_t)
{
    if(!F04TrajPointAllocator)
      F04TrajPointAllocator = new G4Allocator<F04TrajectoryPoint>;
    return (void *) F04TrajPointAllocator->MallocSingle();
}

inline void F04TrajectoryPoint::operator delete(void *aTrajectoryPoint)
{
    F04TrajPointAllocator->FreeSingle(
        (F04TrajectoryPoint *) aTrajectoryPoint);
}

#endif
