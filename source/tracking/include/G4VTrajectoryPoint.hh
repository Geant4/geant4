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
// $Id: G4VTrajectoryPoint.hh,v 1.8 2002-10-16 11:38:37 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4VTrajectoryPoint.hh
//
// class description:
//   This class is the abstract base class which represents the point
//   which consists of a trajectory.
//   It includes information of a the point.
//
// ---------------------------------------------------------------

#ifndef G4VTrajectoryPoint_h
#define G4VTrajectoryPoint_h 1

class G4AttDef;
class G4AttValue;

#include "g4std/vector"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4VTrajectoryPoint
{
 public: // with description

 // Constructor/Destructor
   G4VTrajectoryPoint() {;}
   virtual ~G4VTrajectoryPoint() {;}

 // Operators
   inline int operator==(const G4VTrajectoryPoint& right) const
   { return (this==&right); };

 // Get/Set functions
   virtual const G4ThreeVector GetPosition() const = 0;

 // Get method for a vector of auxiliary points
   virtual const G4std::vector<G4ThreeVector>* GetAuxiliaryPoints() const
   { return 0; }
   // If implemented by a derived class, returns a pointer to a list
   // of auxiliary points, e.g., intermediate points used during the
   // calculation of the step that can be used for drawing a smoother
   // trajectory.  The user must test the validity of this pointer.

 // Get method for HEPRep style attribute definitions
   virtual const G4std::vector<G4AttDef>* GetAttDefs() const
   { return 0; }
   // If implemented by a derived class, returns a pointer to a list
   // of attribute definitions suitable, e.g., for picking.  The user
   // must test the validity of this pointer.

 // Get method for HEPRep style attribute values
   virtual G4std::vector<G4AttValue>* GetAttValues() const
   { return 0; }
   // If implemented by a derived class, returns a pointer to a list
   // of attribute values suitable, e.g., for picking.  The user must
   // test the validity of this pointer and delete the list after use.
};

#endif

