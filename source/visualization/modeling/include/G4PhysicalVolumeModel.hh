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
// $Id: G4PhysicalVolumeModel.hh,v 1.26 2006-03-28 16:50:00 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  31st December 1997.
//
// Class Description:
//
// Model for physical volumes.  It describes a physical volume and its
// daughters to any desired depth.  Note: the "requested depth" is
// specified in the modeling parameters; enum {UNLIMITED = -1}.
//
// For access to base class information, e.g., modeling parameters,
// use GetModelingParameters() inherited from G4VModel.  See Class
// Description of the base class G4VModel.

#ifndef G4PHYSICALVOLUMEMODEL_HH
#define G4PHYSICALVOLUMEMODEL_HH

#include "G4VModel.hh"

#include "G4Transform3D.hh"
#include <iostream>
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSolid;
class G4Material;
class G4VisAttributes;
class G4Polyhedron;

class G4PhysicalVolumeModel: public G4VModel {

public: // With description

  enum {UNLIMITED = -1};

  class G4PhysicalVolumeNodeID {
  public:
    G4PhysicalVolumeNodeID(G4VPhysicalVolume* pPV = 0, G4int iCopyNo = 0):
      fpPV(pPV), fCopyNo(iCopyNo) {}
    G4VPhysicalVolume* GetPhysicalVolume() const {return fpPV;}
    G4int GetCopyNo() const {return fCopyNo;}
    G4bool operator< (const G4PhysicalVolumeNodeID& right) const;
    G4bool operator== (const G4PhysicalVolumeNodeID& right) const;
  private:
    G4VPhysicalVolume* fpPV;
    G4int fCopyNo;
  };
  // Nested class for identifying physical volume nodes.

  G4PhysicalVolumeModel
  (G4VPhysicalVolume*,
   G4int requestedDepth = UNLIMITED,
   const G4Transform3D& modelTransformation = G4Transform3D(),
   const G4ModelingParameters* = 0,
   G4bool useFullExtent = false);

  virtual ~G4PhysicalVolumeModel ();

  void DescribeYourselfTo (G4VGraphicsScene&);
  // The main task of a model is to describe itself to the graphics scene
  // handler (a object which inherits G4VSceneHandler, which inherits
  // G4VGraphicsScene).  It can also provide special information
  // through pointers to working space in the scene handler.  These
  // pointers must be set up (if required by the scene handler) in the
  // scene handler's implementaion of EstablishSpecials
  // (G4PhysicalVolumeModel&) which is called from here.  To do this,
  // the scene handler should call DefinePointersToWorkingSpace - see
  // below.  DecommissionSpecials (G4PhysicalVolumeModel&) is also
  // called from here.  To see how this works, look at the
  // implementation of this function and
  // G4VSceneHandler::Establish/DecommissionSpecials.

  G4String GetCurrentDescription () const;
  // A description which depends on the current state of the model.

  G4String GetCurrentTag () const;
  // A tag which depends on the current state of the model.

  const G4PhysicalVolumeModel* GetG4PhysicalVolumeModel () const
  {return this;}

  G4PhysicalVolumeModel* GetG4PhysicalVolumeModel ()
  {return this;}

  G4VPhysicalVolume* GetTopPhysicalVolume () const {return fpTopPV;}

  G4int GetRequestedDepth () const {return fRequestedDepth;}

  const G4Polyhedron* GetClippingVolume () const
  {return fpClippingPolyhedron;}

  G4int GetCurrentDepth() const {return fCurrentDepth;}
  // Current depth of geom. hierarchy.

  G4VPhysicalVolume* GetCurrentPV() const {return fpCurrentPV;}
  // Current physical volume.

  G4LogicalVolume* GetCurrentLV() const {return fpCurrentLV;}
  // Current logical volume.

  G4Material* GetCurrentMaterial() const {return fpCurrentMaterial;}
  // Current material.

  const std::vector<G4PhysicalVolumeNodeID>& GetDrawnPVPath() const
  {return fDrawnPVPath;}
  // Path of the current drawn (non-culled) volume in terms of drawn
  // (non-culled) ancesters.  It is a vector of physical volume node
  // identifiers corresponding to the geometry hierarchy actually
  // selected, i.e., not culled.  It is guaranteed that volumes are
  // presented to the scene handlers in top-down hierarchy order,
  // i.e., ancesters first, mothers before daughters, so the scene
  // handler can be assured that, if it is building its own scene
  // graph tree, a mother, if any, will have already been encountered
  // and there will already be a node in place on which to hang the
  // current volume.

  void SetRequestedDepth (G4int requestedDepth) {
    fRequestedDepth = requestedDepth;
  }

  void SetClippingPolyhedron (const G4Polyhedron* pClippingPolyhedron) {
    fpClippingPolyhedron = pClippingPolyhedron;
  }

  G4bool Validate (G4bool warn);
  // Validate, but allow internal changes (hence non-const function).

  void DefinePointersToWorkingSpace
  (G4int*              pCurrentDepth = 0,
   G4VPhysicalVolume** ppCurrentPV = 0,
   G4LogicalVolume**   ppCurrentLV = 0,
   G4Material**        ppCurrentMaterial = 0,
   std::vector<G4PhysicalVolumeNodeID>* pDrawnPVPath = 0);
  // For use (optional) by the scene handler if it needs to know about
  // the current information maintained through these pointers.  For
  // description of pointers, see below.
  
  void CurtailDescent() {fCurtailDescent = true;}

protected:

  void VisitGeometryAndGetVisReps (G4VPhysicalVolume*,
				   G4int requestedDepth,
				   const G4Transform3D&,
				   G4VGraphicsScene&);

  void DescribeAndDescend (G4VPhysicalVolume*,
			   G4int requestedDepth,
			   G4LogicalVolume*,
			   G4VSolid*,
			   G4Material*,
			   const G4Transform3D&,
			   G4VGraphicsScene&);

  virtual void DescribeSolid (const G4Transform3D& theAT,
			      G4VSolid* pSol,
			      const G4VisAttributes* pVisAttribs,
			      G4VGraphicsScene& sceneHandler);

  void CalculateExtent ();

  /////////////////////////////////////////////////////////
  // Data members...

  G4VPhysicalVolume* fpTopPV;        // The physical volume.
  G4String           fTopPVName;     // ...of the physical volume.
  G4int              fTopPVCopyNo;   // ...of the physical volume.
  G4int              fRequestedDepth;
                     // Requested depth of geom. hierarchy search.
  G4bool             fUseFullExtent; // ...if requested.
  G4int              fCurrentDepth;  // Current depth of geom. hierarchy.
  G4VPhysicalVolume* fpCurrentPV;    // Current physical volume.
  G4LogicalVolume*   fpCurrentLV;    // Current logical volume.
  G4Material*    fpCurrentMaterial;  // Current material.
  std::vector<G4PhysicalVolumeNodeID> fDrawnPVPath;
  G4bool             fCurtailDescent;  // Can be set to curtail descent.
  const G4Polyhedron*fpClippingPolyhedron;

  ////////////////////////////////////////////////////////////
  // Pointers to working space in scene handler, if required.

  G4int*              fpCurrentDepth;  // Current depth of geom. hierarchy.
  G4VPhysicalVolume** fppCurrentPV;    // Current physical volume.
  G4LogicalVolume**   fppCurrentLV;    // Current logical volume.
  G4Material**    fppCurrentMaterial;  // Current material.
  std::vector<G4PhysicalVolumeNodeID>* fpDrawnPVPath;

private:

  // Private copy constructor and assigment operator - copying and
  // assignment not allowed.  Keeps CodeWizard happy.
  G4PhysicalVolumeModel (const G4PhysicalVolumeModel&);
  G4PhysicalVolumeModel& operator = (const G4PhysicalVolumeModel&);
};

std::ostream& operator<<
  (std::ostream& os, const G4PhysicalVolumeModel::G4PhysicalVolumeNodeID);

#endif
