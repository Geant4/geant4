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
//
// G4PhysicalVolumeModel assumes the modeling parameters have been set
// up with meaningful information - default vis attributes and culling
// policy in particular.
//
// The volumes are unpacked and sent to the scene handler.  They are,
// in effect, "touchables" as defined by G4TouchableHistory.  A
// touchable is defined not only by its physical volume name and copy
// number but also by its position in the geometry hierarchy.
//
// It is guaranteed that touchables are presented to the scene handler
// in top-down hierarchy order, i.e., ancestors first, mothers before
// daughters, so the scene handler can be assured that, if it is
// building its own scene graph tree, a mother, if any, will have
// already been encountered and there will already be a node in place
// on which to hang the current volume.  But be aware that the
// visibility and culling policy might mean that some touchables are
// not passed to the scene handler so the drawn tree might have
// missing layers.  GetFullPVPath allows you to know the full
// hierarchy.

#ifndef G4PHYSICALVOLUMEMODEL_HH
#define G4PHYSICALVOLUMEMODEL_HH

#include "G4VModel.hh"
#include "G4ModelingParameters.hh"

#include "G4VTouchable.hh"
#include "G4Transform3D.hh"
#include "G4Plane3D.hh"
#include <iostream>
#include <vector>
#include <map>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSolid;
class G4Material;
class G4VisAttributes;
class G4AttDef;
class G4AttValue;

class G4PhysicalVolumeModel: public G4VModel {

public: // With description

  enum {UNLIMITED = -1};

  enum ClippingMode {subtraction, intersection};

  // Nested class for identifying physical volume nodes.
  class G4PhysicalVolumeNodeID {
  public:
    G4PhysicalVolumeNodeID
    (G4VPhysicalVolume* pPV = 0,
     G4int iCopyNo = 0,
     G4int depth = 0,
     const G4Transform3D& transform = G4Transform3D(),
     G4bool drawn = true):
      fpPV(pPV),
      fCopyNo(iCopyNo),
      fNonCulledDepth(depth),
      fTransform(transform),
      fDrawn(drawn) {}
    G4VPhysicalVolume* GetPhysicalVolume() const {return fpPV;}
    G4int  GetCopyNo()                     const {return fCopyNo;}
    G4int  GetNonCulledDepth()             const {return fNonCulledDepth;}
    const  G4Transform3D& GetTransform()   const {return fTransform;}
    G4bool GetDrawn()                      const {return fDrawn;}
    void SetPhysicalVolume(G4VPhysicalVolume* v)   {fpPV = v;}
    void SetCopyNo        (G4int n)                {fCopyNo = n;}
    void SetNonCulledDepth(G4int d)                {fNonCulledDepth = d;}
    void SetTransform     (const G4Transform3D& t) {fTransform = t;}
    void SetDrawn         (G4bool b)               {fDrawn = b;}
    G4bool operator<  (const G4PhysicalVolumeNodeID& right) const;
    G4bool operator!= (const G4PhysicalVolumeNodeID& right) const;
    G4bool operator== (const G4PhysicalVolumeNodeID& right) const {
      return !operator!= (right);
    }
  private:
    G4VPhysicalVolume* fpPV;
    G4int fCopyNo;
    G4int fNonCulledDepth;
    G4Transform3D fTransform;
    G4bool fDrawn;
  };

  // Nested class for handling nested parameterisations.
  class G4PhysicalVolumeModelTouchable: public G4VTouchable {
  public:
    G4PhysicalVolumeModelTouchable
    (const std::vector<G4PhysicalVolumeNodeID>& fullPVPath);
    const G4ThreeVector& GetTranslation(G4int depth) const;
    const G4RotationMatrix* GetRotation(G4int depth) const;
    G4VPhysicalVolume* GetVolume(G4int depth) const;
    G4VSolid* GetSolid(G4int depth) const;
    G4int GetReplicaNumber(G4int depth) const;
    G4int GetHistoryDepth() const {return G4int(fFullPVPath.size());}
  private:
    const std::vector<G4PhysicalVolumeNodeID>& fFullPVPath;
  };

  // Nested struct for encapsulating touchable properties
  struct TouchableProperties {
    TouchableProperties(): fpTouchablePV(nullptr), fCopyNo(0) {}
    G4ModelingParameters::PVNameCopyNoPath fTouchablePath;
    G4VPhysicalVolume*                     fpTouchablePV;
    G4int                                  fCopyNo;
    G4Transform3D                          fTouchableGlobalTransform;
    std::vector<G4PhysicalVolumeNodeID>    fTouchableBaseFullPVPath;
    std::vector<G4PhysicalVolumeNodeID>    fTouchableFullPVPath;
  };

  G4PhysicalVolumeModel
  (G4VPhysicalVolume* = 0,
   G4int requestedDepth = UNLIMITED,
   const G4Transform3D& modelTransformation = G4Transform3D(),
   const G4ModelingParameters* = 0,
   G4bool useFullExtent = false,
   const std::vector<G4PhysicalVolumeNodeID>& baseFullPVPath =
   std::vector<G4PhysicalVolumeNodeID>());

  virtual ~G4PhysicalVolumeModel ();

  void DescribeYourselfTo (G4VGraphicsScene&);
  // The main task of a model is to describe itself to the graphics scene
  // handler (a object which inherits G4VSceneHandler, which inherits
  // G4VGraphicsScene).

  G4String GetCurrentDescription () const;
  // A description which depends on the current state of the model.

  G4String GetCurrentTag () const;
  // A tag which depends on the current state of the model.

  G4VPhysicalVolume* GetTopPhysicalVolume () const {return fpTopPV;}

  G4int GetRequestedDepth () const {return fRequestedDepth;}

  const G4VSolid* GetClippingSolid () const
  {return fpClippingSolid;}

  G4int GetCurrentDepth() const {return fCurrentDepth;}
  // Current depth in geom. hierarchy.

  const G4Transform3D& GetTransformation() const {return fTransform;}
  // Initial transformation

  G4VPhysicalVolume* GetCurrentPV() const {return fpCurrentPV;}
  // Current physical volume.

  G4int GetCurrentPVCopyNo() const {return fCurrentPVCopyNo;}
  // Current copy number.

  G4LogicalVolume* GetCurrentLV() const {return fpCurrentLV;}
  // Current logical volume.

  G4Material* GetCurrentMaterial() const {return fpCurrentMaterial;}
  // Current material.

  const G4Transform3D& GetCurrentTransform() const {return fCurrentTransform;}
  // Current transform.

  const std::vector<G4PhysicalVolumeNodeID>& GetBaseFullPVPath() const
  {return fBaseFullPVPath;}
  // Base for this G4PhysicalVolumeModel. It can be empty, which would be the
  // case for the world volume. But the user may create a G4PhysicalVolumeModel
  // for a volume anywhere in the tree, in which case the user must provide
  // the transform (in the constructor) and the base path (SetBaseFullPVPath).

  const std::vector<G4PhysicalVolumeNodeID>& GetFullPVPath() const
  {return fFullPVPath;}
  // Vector of physical volume node identifiers for the current
  // touchable.  It is its path in the geometry hierarchy, similar to
  // the concept of "touchable history" available from the navigator
  // during tracking.

  const std::vector<G4PhysicalVolumeNodeID>& GetDrawnPVPath() const
  {return fDrawnPVPath;}
  // Path of the current drawn (non-culled) touchable in terms of
  // drawn (non-culled) ancestors.  It is a vector of physical volume
  // node identifiers corresponding to the geometry hierarchy actually
  // selected, i.e., with "culled" volumes NOT included.

  static G4ModelingParameters::PVNameCopyNoPath GetPVNameCopyNoPath
  (const std::vector<G4PhysicalVolumeNodeID>&);
  // Converts

  const std::map<G4String,G4AttDef>* GetAttDefs() const;
  // Attribute definitions for current solid.

  std::vector<G4AttValue>* CreateCurrentAttValues() const;
  // Attribute values for current solid.  Each must refer to an
  // attribute definition in the above map; its name is the key.  The
  // user must test the validity of this pointer (it must be non-zero
  // and conform to the G4AttDefs, which may be checked with
  // G4AttCheck) and delete the list after use.  See
  // G4XXXStoredSceneHandler::PreAddSolid for how to access and
  // G4VTrajectory::ShowTrajectory for an example of the use of
  // G4Atts.

  void SetRequestedDepth (G4int requestedDepth) {
    fRequestedDepth = requestedDepth;
  }

  void SetClippingSolid (G4VSolid* pClippingSolid) {
    fpClippingSolid = pClippingSolid;
  }

  void SetClippingMode (ClippingMode mode) {
    fClippingMode = mode;
  }

  G4bool Validate (G4bool warn);
  // Validate, but allow internal changes (hence non-const function).

  void Abort () const {fAbort = true;}
  // Abort all further traversing.

  void CurtailDescent() const {fCurtailDescent = true;}
  // Curtail descent of current branch.

  void CalculateExtent ();

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

  /////////////////////////////////////////////////////////
  // Data members...

  G4VPhysicalVolume* fpTopPV;        // The physical volume.
  G4String           fTopPVName;     // ...of the physical volume.
  G4int              fTopPVCopyNo;   // ...of the physical volume.
  G4int              fRequestedDepth;
                     // Requested depth of geom. hierarchy search.
  G4bool             fUseFullExtent; // ...if requested.
  G4Transform3D      fTransform;     // Initial transform
  G4int              fCurrentDepth;  // Current depth of geom. hierarchy.
  G4VPhysicalVolume* fpCurrentPV;    // Current physical volume.
  G4int              fCurrentPVCopyNo; // Current copy number.
  G4LogicalVolume*   fpCurrentLV;    // Current logical volume.
  G4Material*        fpCurrentMaterial;  // Current material.
  G4Transform3D      fCurrentTransform;  // Current transform.
  std::vector<G4PhysicalVolumeNodeID> fBaseFullPVPath; // Base. May be empty.
  std::vector<G4PhysicalVolumeNodeID> fFullPVPath;     // Starts from base.
  std::vector<G4PhysicalVolumeNodeID> fDrawnPVPath;    // Omits culled volumes.
  mutable G4bool     fAbort;         // Abort all further traversing.
  mutable G4bool     fCurtailDescent;// Can be set to curtail descent.
  G4VSolid*          fpClippingSolid;
  ClippingMode       fClippingMode;

private:

  // Private copy constructor and assigment operator - copying and
  // assignment not allowed.  Keeps CodeWizard happy.
  G4PhysicalVolumeModel (const G4PhysicalVolumeModel&);
  G4PhysicalVolumeModel& operator = (const G4PhysicalVolumeModel&);
};

std::ostream& operator<<
(std::ostream& os, const G4PhysicalVolumeModel::G4PhysicalVolumeNodeID&);

std::ostream& operator<<
(std::ostream& os, const std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>&);

#endif
