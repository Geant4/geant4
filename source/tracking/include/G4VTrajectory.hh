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
// $Id: G4VTrajectory.hh,v 1.11 2002-11-08 18:17:04 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4VTrajectory.hh
//
// class description:
//   This class is the abstract base class which represents a trajectory of 
//   a particle tracked.
//   Its concrete class includes information of 
//     1) List of trajectory points which compose the trajectory,
//     2) static information of particle which generated the 
//        trajectory, 
//     3) trackID and parent particle ID of the trajectory,
//
// ---------------------------------------------------------------

#ifndef G4VTrajectory_h
#define G4VTrajectory_h 1

#include "globals.hh"
#include "g4std/vector"
#include "g4std/map"
#include "G4ThreeVector.hh"

class G4Step;
class G4VTrajectoryPoint;
class G4AttDef;
class G4AttValue;

class G4VTrajectory
{
 public: // with description

// Constructor/Destrcutor

   G4VTrajectory() {;}
   virtual ~G4VTrajectory() {;}

// Operators
   inline int operator == (const G4VTrajectory& right) const
   {return (this==&right);} 

// Get/Set functions 
   virtual G4int GetTrackID() const = 0;
   virtual G4int GetParentID() const = 0;
   virtual G4String GetParticleName() const = 0;
   virtual G4double GetCharge() const = 0;
   // Charge is that of G4DynamicParticle
   virtual G4int GetPDGEncoding() const = 0;
   // Zero will be returned if the particle does not have PDG code.
   virtual G4ThreeVector GetInitialMomentum() const = 0;
   // Momentum at the origin of the track in global coordinate system.

// Other member functions
   virtual int GetPointEntries() const = 0;
   // Returns the number of trajectory points
   virtual G4VTrajectoryPoint* GetPoint(G4int i) const = 0;
   // Returns i-th trajectory point.
   virtual void ShowTrajectory(G4std::ostream& os=G4cout) const;
   // Convert attributes in trajectory (and trajectory point if
   // needed) to ostream.  A default implementation in this base class
   // may be used or may be overridden in the concrete class.  Note:
   // the user needs to follow with new-line or end-of-string,
   // depending on the nature of os.
   virtual void DrawTrajectory(G4int i_mode=0) const;
   // Draw the trajectory.  A default implementation in this base
   // class may be used or may be overridden in the concrete class.
   virtual const G4std::map<G4String,G4AttDef>* GetAttDefs() const
   { return 0; }
   // If implemented by a derived class, returns a pointer to a map of
   // attribute definitions for the attribute values below.  The user
   // must test the validity of this pointer.
   virtual G4std::vector<G4AttValue>* CreateAttValues() const
   { return 0; }
   // If implemented by a derived class, returns a pointer to a list
   // of attribute values suitable, e.g., for picking.  Each must
   // refer to an attribute definition in the above map; its name is
   // the key.  The user must test the validity of this pointer and
   // delete the list after use.

 public:
   // Following methods MUST be invoked exclusively by G4TrackingManager
   virtual void AppendStep(const G4Step* aStep) = 0;
   virtual void MergeTrajectory(G4VTrajectory* secondTrajectory) = 0;

};

#endif










