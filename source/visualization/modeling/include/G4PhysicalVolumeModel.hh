// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicalVolumeModel.hh,v 1.6 1999-11-25 14:18:49 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  31st December 1997.
//
// Class Description:
//
// Model for physical volumes.  It describes a physical volume and its
// daughters to any desired depth.  Note: the "sought depth" is
// specified in the modeling parameters.
//
// For access to base class information, e.g., modeling parameters,
// use GetModelingParameters() inherited from G4VModel.  See Class
// Description of the base class G4VModel.

#ifndef G4PHYSICALVOLUMEMODEL_HH
#define G4PHYSICALVOLUMEMODEL_HH

#include "G4VModel.hh"

#include "G4Transform3D.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSolid;
class G4Material;
class G4VisAttributes;

class G4PhysicalVolumeModel: public G4VModel {

public: // With description

  enum {UNLIMITED = -1};

  G4PhysicalVolumeModel
  (G4VPhysicalVolume*,
   G4int soughtDepth = UNLIMITED,
   const G4Transform3D& modelTransformation = G4Transform3D::Identity,
   const G4ModelingParameters* = 0,
   G4bool useFullExtent = false);

  virtual ~G4PhysicalVolumeModel ();

  void DescribeYourselfTo (G4VGraphicsScene&);
  // The main task of a model is to describe itself to the scene
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

  G4bool Validate ();
  // Validate, but allow internal changes (hence non-const function).

  void DefinePointersToWorkingSpace (G4int*              pCurrentDepth,
				     G4VPhysicalVolume** ppCurrentPV,
				     G4LogicalVolume**   ppCurrentLV);
  // For use (optional) by the scene handler if it needs to know about
  // the current information maintained through these pointers.

protected:

  void VisitGeometryAndGetVisReps (G4VPhysicalVolume*,
				   G4int soughtDepth,
				   const G4Transform3D&,
				   G4VGraphicsScene&);

  void DescribeAndDescend (G4VPhysicalVolume*,
			   G4int soughtDepth,
			   G4LogicalVolume*,
			   G4VSolid*,
			   const G4Material*,
			   const G4Transform3D&,
			   G4VGraphicsScene&);

  void DescribeSolids (const G4Transform3D& theAT,
		       G4VSolid* pSol,
		       const G4VisAttributes* pVisAttribs,
		       G4VGraphicsScene& sceneHandler);

  G4bool IsThisCulled     (const G4LogicalVolume*,
			   const G4Material*);
  G4bool IsDaughterCulled (const G4LogicalVolume* pMotherLV);

  /////////////////////////////////////////////////////////
  // Data members...

  G4VPhysicalVolume* fpTopPV;        // The physical volume.
  G4String           fTopPVName;     // ...of the physical volume.
  G4int              fTopPVCopyNo;   // ...of the physical volume.
  G4int              fSoughtDepth;   // Sought depth of geom. hierarchy search.
  G4int              fCurrentDepth;  // Current depth of geom. hierarchy.
  G4VPhysicalVolume* fpCurrentPV;    // Current physical volume.
  G4LogicalVolume*   fpCurrentLV;    // Current logical volume.

  ////////////////////////////////////////////////////////////
  // Pointers to working space in scene, if required.

  G4int*              fpCurrentDepth;  // Current depth of geom. hierarchy.
  G4VPhysicalVolume** fppCurrentPV;    // Current physical volume.
  G4LogicalVolume**   fppCurrentLV;    // Current logical volume.

};

#endif
