// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ViewParameters.hh,v 1.1 1999-01-07 16:15:19 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  19th July 1996
// View parameters and options.

#ifndef G4VIEWPARAMETERS_HH
#define G4VIEWPARAMETERS_HH

#include <rw/tvordvec.h>
#include "G4Vector3D.hh"
#include "G4Point3D.hh"
#include "G4Plane3D.hh"
#include "G4VisAttributes.hh"
#include "G4VMarker.hh"

typedef RWTValOrderedVector<G4Plane3D> G4Planes;

//. % THE STANDARD VIEW AND ALL THAT.
//.
//. In GEANT4 visualization, we have the concept of a ``Standard View''.
//. This is the view when the complete set of objects being viewed is
//. comfortably in view from any viewpoint.  It is defined by the
//. ``Bounding Sphere'' of ``visible'' objects when initially registered
//. in the scene, and the View Parameters.
//.
//. There is also the ``Standard Target Point'', which is the centre
//. of the Bounding Sphere (note that these belong to the scene, and are
//. stored in the Scene Data).  The ``Current Target Point'', initially
//. set to the Standard Target Point, is incremented by the
//. ``dolly'' and ``zoom'' commands, and can be reset to the
//. Standard Target Point with the {\tt /vis~/camera/reset} command.
//.
//. Also, the ``Standard Camera Position'' is the ``Standard Camera
//. Distance'' along
//. the Viewpoint Direction vector from the Standard Target Point.
//. The Standard Camera Distance is the radius of the
//. Bounding Sphere divided by {\tt fFieldHalfAngle}.  It is not stored
//. explicitly because of the singularity at {\tt fFieldHalfAngle} = 0,
//. which implies parallel projection.
//.
//. Similarly, the ``Current Camera Position'' is the ``Current
//. Camera Distance''
//. along the Viewpoint Direction vector from the Current Target Point.
//. The Current Camera Distance is given by the formulae below, but
//. note that it can be negative, meaning that the camera has
//. moved {\em beyond} the Current Target Point, which is conceptually
//. possible, but which might give some problems when setting up
//. the view matrix --- see, for example, {\tt G4OpenGLView::SetView ()}.
//.
//. All views are expected to keep the ``Up Vector'' vertical.
//.
//. Finally, the view is magnified by the ``Zoom Factor'' which is
//. reset to 1 by the \linebreak
//. {\tt /vis~/camera/reset} command.
//.
//. Note: the camera movements pan, dolly and zoom, are currently coded in
//. \linebreak {\tt G4VisManMessenger.cc}.
//.
//. The algorithms for calculating various useful quantities from the
//. View Parameters are encapsulated in the functions described in
//. Appendix \ref{ap:useful}.

class G4ViewParameters {

  friend ostream& operator << (ostream& os, const G4ViewParameters& v);
  friend G4bool   operator != (const G4ViewParameters& v1,
			       const G4ViewParameters& v2);
public:

  enum DrawingStyle {
    wireframe,  // Draw edges    - no hidden line removal.
    hlr,        // Draw edges    - hidden lines removed.
    hsr,        // Draw surfaces - hidden surfaces removed.
    hlhsr       // Draw surfaces and edges - hidden removed.
  };

  enum RepStyle {
    polyhedron, // Use G4Polyhedron.
    nurbs       // Use G4NURBS.
  };

  G4ViewParameters ();
  ~G4ViewParameters ();

  // Note: uses default assignment operator and copy constructor.

  // Get and Is functions.
        DrawingStyle     GetDrawingStyle         () const;
        RepStyle         GetRepStyle             () const;
        G4bool           IsCulling               () const;
        G4bool           IsCullingInvisible      () const;
        G4bool           IsDensityCulling        () const;
        G4double         GetVisibleDensity       () const;
        G4bool           IsCullingCovered        () const;
        G4bool           IsSection               () const;
  const G4Plane3D&       GetSectionPlane         () const;
        G4bool           IsCutaway               () const;
  const G4Planes&        GetCutawayPlanes        () const;
        G4bool           IsExplode               () const;
        G4double         GetExplodeFactor        () const;
        G4int            GetNoOfSides            () const;
  const G4Vector3D&      GetViewpointDirection   () const;
  const G4Vector3D&      GetUpVector             () const;
        G4double         GetFieldHalfAngle       () const;
        G4double         GetZoomFactor           () const;
  const G4Point3D&       GetCurrentTargetPoint   () const;
        G4double         GetDolly                () const;
        G4bool           GetLightsMoveWithCamera () const;
  const G4Vector3D&      GetLightpointDirection  () const;  // Relative...
        G4Vector3D&      GetActualLightpointDirection  ();  // Actual...
  // ... depending on GetLightsMoveWithCamera.
        G4bool           IsViewGeom              () const;
        G4bool           IsViewHits              () const;
        G4bool           IsViewDigis             () const;
  const G4VisAttributes* GetDefaultVisAttributes () const;
  const G4VisAttributes* GetDefaultTextVisAttributes () const;
  const G4VMarker&       GetDefaultMarker        () const;
        G4double         GetGlobalMarkerScale    () const;
        G4bool           IsMarkerNotHidden       () const;
        G4int            GetWindowSizeHintX      () const;
        G4int            GetWindowSizeHintY      () const;

  // Here Follow functions to evaluate the above algorithms as a
  // function of the radius of the Bounding Sphere of the object being
  // viewed.  Call them in the order given - for efficiency, later
  // functions depend on the results of earlier ones (Store the
  // results of earlier functions in your own temporary variables -
  // see, for example, G4OpenGLView::SetView ().)
  G4double GetCameraDistance  (G4double radius) const;
  G4double GetNearDistance    (G4double cameraDistance, G4double radius) const;
  G4double GetFarDistance     (G4double cameraDistance,
			       G4double nearDistance, G4double radius) const;
  G4double GetFrontHalfHeight (G4double nearDistance, G4double radius) const;

  // Set, Add, Multiply, Increment, Unset and Clear functions.
  void SetDrawingStyle         (G4ViewParameters::DrawingStyle style);
  void SetRepStyle             (G4ViewParameters::RepStyle style);
  void SetCulling              (G4bool);
  void SetCullingInvisible     (G4bool);
  void SetDensityCulling       (G4bool);
  void SetVisibleDensity       (G4double visibleDensity);
  void SetCullingCovered       (G4bool);
  void SetSectionPlane         (const G4Plane3D& sectionPlane);
  void UnsetSectionPlane       ();
  void AddCutawayPlane         (const G4Plane3D& cutawayPlane);
  void ClearCutawayPlanes      ();
  void SetExplodeFactor        (G4double explodeFactor);
  void UnsetExplodeFactor      ();
  void SetNoOfSides            (G4int nSides);
  void SetViewpointDirection   (const G4Vector3D& viewpointDirection);
  // Prefer the following to get lightpoint direction right too.
  void SetViewAndLights        (const G4Vector3D& viewpointDirection);
  // Also sets lightpoint direction according to G4bool fLightsMoveWithCamera.
  void SetUpVector             (const G4Vector3D& upVector);
  void SetFieldHalfAngle       (G4double fieldHalfAngle);
  void SetZoomFactor           (G4double zoomFactor);
  void MultiplyZoomFactor      (G4double zoomFactorMultiplier);
  void SetCurrentTargetPoint   (const G4Point3D& currentTargetPoint);
  void SetDolly                (G4double dolly);
  void IncrementDolly          (G4double dollyIncrement);
  void SetLightpointDirection  (const G4Vector3D& lightpointDirection);
  void SetLightsMoveWithCamera (G4bool moves);
  void SetViewGeom             ();
  void UnsetViewGeom           ();
  void SetViewHits             ();
  void UnsetViewHits           ();
  void SetViewDigis            ();
  void UnsetViewDigis          ();
  void SetDefaultVisAttributes (const G4VisAttributes&);
  void SetDefaultTextVisAttributes (const G4VisAttributes&);
  void SetDefaultMarker        (const G4VMarker& defaultMarker);
  void SetGlobalMarkerScale    (G4double globalMarkerScale);
  void SetMarkerHidden         ();
  void SetMarkerNotHidden      ();
  void SetWindowSizeHint       (G4int xHint, G4int yHint);

  void Pan                     (G4double right, G4double up);
  // Note: the result of above "pan"operation is always an Increment
  // relative to the current target point and also depends on the
  // veiwpoint and up directions.

  void PrintDifferences (const G4ViewParameters& v) const;

private:

  DrawingStyle fDrawingStyle;    // Drawing style.
  RepStyle     fRepStyle;        // Representation style.
  G4bool       fCulling;         // Culling requested.
  G4bool       fCullInvisible;   // Cull (don't Draw) invisible objects.
  G4bool       fDensityCulling;  // Density culling requested.  If so...
  G4double     fVisibleDensity;  // ...density lower than this not drawn.
  G4bool       fCullCovered;     // Cull daughters covered by opaque mothers.
  G4bool       fSection;         // Section drawing requested (DCUT in GEANT3).
  G4Plane3D    fSectionPlane;    // Cut plane for section drawing (DCUT).
  G4bool       fCutaway;         // Cutaway flag.
  G4Planes     fCutawayPlanes;   // Set of planes used for cutaway.
  G4bool       fExplode;         // Explode flag.
  G4double     fExplodeFactor;
  G4int        fNoOfSides;       // ...if polygon approximates circle.
  G4Vector3D   fViewpointDirection;
  G4Vector3D   fUpVector;        // Up vector.  (Warning: MUST NOT be parallel
                                 // to fViewpointDirection!)
  G4double     fFieldHalfAngle;  // Radius / camara distance, 0 for parallel.
  G4double     fZoomFactor;      // Magnification relative to Standard View.
  G4Point3D    fCurrentTargetPoint;
  G4double     fDolly;           // Distance towards current target point.
  G4bool       fLightsMoveWithCamera;
  G4Vector3D   fRelativeLightpointDirection;
  // i.e., rel. to object or camera accoding to G4bool LightsMoveWithCamera.
  G4Vector3D   fActualLightpointDirection;
  G4bool       fViewGeom;        // View geometry objects.
  G4bool       fViewHits;        // View hits, if any.
  G4bool       fViewDigis;       // View digis, if any.
  G4VisAttributes fDefaultVisAttributes;
  G4VisAttributes fDefaultTextVisAttributes;
  G4VMarker    fDefaultMarker;
  G4double     fGlobalMarkerScale;
  G4bool       fMarkerNotHidden;
  // True if transients are to be drawn and not hidden by
  // hidden-line-hidden-surface removal algorithms, e.g., z-buffer
  // testing; false if they are to be hidden-line-hidden-surface
  // removed.
  G4int        fWindowSizeHintX; // Size hints for pixel-based window systems.
  G4int        fWindowSizeHintY;
};

#include "G4ViewParameters.icc"

#endif

