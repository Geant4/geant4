// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VSceneHandler.hh,v 1.7 1999-11-11 15:38:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  19th July 1996.
// Abstract interface class for graphics scene handlers.
// Inherits from G4VGraphicsScene, in the graphics_reps component, which is
// a minimal abstract interface for the GEANT4 kernel.

#ifndef G4VSCENEHANDLER_HH
#define G4VSCENEHANDLER_HH

#include "globals.hh"

#include "G4VGraphicsScene.hh"
#include "G4Scene.hh"
#include "G4ViewerList.hh"
#include "G4ViewParameters.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"
#include "G4Torus.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4Transform3D.hh"
#include "G4VModel.hh"

class G4VGraphicsSystem;
class G4VViewer;
class G4VSolid;
class G4VisAttributes;
class G4Colour;
class G4Visible;
class G4VSolid;
class G4ModelingParameters;
class G4LogicalVolume;

class G4VSceneHandler: public G4VGraphicsScene {

  friend ostream& operator << (ostream& os, const G4VSceneHandler& s);

public:

  enum MarkerSizeType {world, screen};

  G4VSceneHandler (G4VGraphicsSystem& system, G4int id, const G4String& name = "");

  virtual ~G4VSceneHandler ();

  // For G4RWTPtrOrderedVector...
  G4bool operator == (const G4VSceneHandler& scene) const;

  //////////////////////////////////////////////////////////////
  // Functions for adding raw GEANT4 objects, if the graphics system
  // can can understand them (optional on the part of the graphics
  // system).  If your graphics system is sophisticated enough to
  // handle a particular solid shape as a primitive, in your derived
  // class write a function to override one or more of the following.
  // See the implementation of G4VSceneHandler::AddThis (const G4Box& box) for
  // more suggestions.  If not, please implement the base class
  // invocation.
  virtual void AddThis (const G4Box&);
  virtual void AddThis (const G4Cons&);
  virtual void AddThis (const G4Tubs&);
  virtual void AddThis (const G4Trd&);
  virtual void AddThis (const G4Trap&);
  virtual void AddThis (const G4Sphere&);
  virtual void AddThis (const G4Para&);
  virtual void AddThis (const G4Torus&);
  virtual void AddThis (const G4Polycone&);
  virtual void AddThis (const G4Polyhedra&);
  virtual void AddThis (const G4VSolid&);  // For solids not above.

  ///////////////////////////////////////////////////////////////
  // Other inherited functions.

  virtual void EstablishSpecials (G4PhysicalVolumeModel&);
  // Used to establish any special relationships between scene and this
  // particular type of model - non-pure, i.e., no requirement to
  // implement.  See G4PhysicalVolumeModel.hh for details.

  virtual void PreAddThis (const G4Transform3D& objectTransformation,
			   const G4VisAttributes& visAttribs);
  // objectTransformation is the transformation in the world
  // coordinate system of the object about to be added, and visAttribs
  // is its visualization attributes.
  // IMPORTANT: invoke this from your polymorphic versions, e.g.:
  // void MyXXXSceneHandler::PreAddThis
  //  (const G4Transform3D& objectTransformation,
  //   const G4VisAttributes& visAttribs) {
  //   G4VSceneHandler::PreAddThis (objectTransformation, visAttribs);
  //   ...
  // }

  virtual void PostAddThis ();
  // IMPORTANT: invoke this from your polymorphic versions, e.g.:
  // void MyXXXSceneHandler::PostAddThis () {
  //   ...
  //   G4VSceneHandler::PostAddThis (objectTransformation, visAttribs);
  // }

  //////////////////////////////////////////////////////////////
  // Functions for adding primitives.

  virtual void BeginModeling ();
  // IMPORTANT: invoke this from your polymorphic versions, e.g.:
  // void MyXXXScene::BeginModeling () {
  //   G4VSceneHandler::BeginModeling ();
  //   ...
  // }

  virtual void EndModeling ();
  // IMPORTANT: invoke this from your polymorphic versions, e.g.:
  // void MyXXXScene::EndModeling () {
  //   ...
  //   G4VSceneHandler::EndModeling ();
  // }

  virtual void BeginPrimitives (const G4Transform3D& objectTransformation);
  // IMPORTANT: invoke this from your polymorphic versions, e.g.:
  // void MyXXXSceneHandler::BeginPrimitives
  // (const G4Transform3D& objectTransformation) {
  //   G4VSceneHandler::BeginPrimitives (objectTransformation);
  //   ...
  // }

  virtual void EndPrimitives ();
  // IMPORTANT: invoke this from your polymorphic versions, e.g.:
  // void MyXXXSceneHandler::EndPrimitives () {
  //   ...
  //   G4VSceneHandler::EndPrimitives ();
  // }

  virtual void AddPrimitive (const G4Polyline&)   = 0;
  virtual void AddPrimitive (const G4Text&)       = 0;
  virtual void AddPrimitive (const G4Circle&)     = 0;      
  virtual void AddPrimitive (const G4Square&)     = 0;      
  virtual void AddPrimitive (const G4Polymarker&);  // Implemented in terms
  // of the above but can be over-ridden.
  virtual void AddPrimitive (const G4Polyhedron&) = 0;  
  virtual void AddPrimitive (const G4NURBS&)      = 0;       

  //////////////////////////////////////////////////////////////
  // Access functions.
  const G4String&     GetName           () const;
  void                SetName           (const G4String&);
  G4int               GetSceneId        () const;
  G4int               GetViewCount      () const;
  G4VGraphicsSystem*  GetGraphicsSystem () const;
  G4Scene*            GetScene          () const;
  const G4ViewerList& GetViewerList     () const;
  const G4VModel*     GetModel          () const;
  G4VViewer*          GetCurrentViewer  () const;
  void          SetCurrentViewer (G4VViewer*);
  void          SetScene         (G4Scene*);
  G4ViewerList& SetViewerList    ();  // Non-const so you can change.
  void          SetModel         (const G4VModel*);

  //////////////////////////////////////////////////////////////
  // Public utility functions.

  const G4Colour& GetColour (const G4Visible&);
  const G4Colour& GetColor  (const G4Visible&);
  // The above return colour of G4Visible object, or default global colour.

  const G4Colour& GetTextColour (const G4Text&);
  const G4Colour& GetTextColor  (const G4Text&);
  // The above return colour of G4Text object, or default text colour.

  G4ViewParameters::DrawingStyle GetDrawingStyle (const G4Visible&);
  // Returns drawing style of G4Visible object, i.e., the global
  // default value unless ovewrridden by forced attributes.  And the
  // attributes themselves can be specified or the global default is
  // used - see next function.  So use this function is you have a
  // G4Visible, or the next if all you have is a G4VisAttributes
  // pointer.

  G4ViewParameters::DrawingStyle GetDrawingStyle (const G4VisAttributes*);
  // Returns global default drawing style unless forced attributes are set.
  // So always use this to get the applicable drawing style.

  G4double GetMarkerSize (const G4VMarker&, MarkerSizeType&);
  // Returns applicable marker size (diameter) and type (in second
  // argument).  Uses global default marker if marker sizes are not
  // set.

  G4double GetMarkerDiameter (const G4VMarker&, MarkerSizeType&);
  // Alias for GetMarkerSize.

  G4double GetMarkerRadius (const G4VMarker&, MarkerSizeType&);
  // GetMarkerSize / 2.

  G4ModelingParameters* CreateModelingParameters ();
  // Only the scene and view know what the Modeling Parameters should
  // be.  For historical reasons, the GEANT4 Visualization Environment
  // maintains its own Scene Data and View Parameters, which must be
  // converted, when needed, to Modeling Parameters.

  //////////////////////////////////////////////////////////////
  // Administration functions.

  G4int IncrementViewCount ();

  virtual void ClearStore ();
  // Clears graphics database (display lists) if any.  This base class
  // implements some common functionality so...
  // IMPORTANT: invoke this from your polymorphic versions, e.g.:
  // void MyXXXSceneHandler::ClearStore () {
  //   G4VSceneHandler::ClearStore ();
  //   ...
  // }

  virtual void ClearTransientStore ();
  // Clears transient part of graphics database (display lists) if any.
  // This base class implements some common functionality so...
  // IMPORTANT: invoke this from your polymorphic versions, e.g.:
  // void MyXXXSceneHandler::ClearTransientStore () {
  //   G4VSceneHandler::ClearTransientStore ();
  //   ...
  // }

  void AddViewerToList      (G4VViewer* pView);  // Add view to view List.
  void RemoveViewerFromList (G4VViewer* pView);  // Remove view from view List.

protected:

  //////////////////////////////////////////////////////////////
  // Default routine used by default AddThis ().

  virtual void RequestPrimitives (const G4VSolid& solid);

  //////////////////////////////////////////////////////////////
  // Data members

  G4VGraphicsSystem&     fSystem;      // Graphics system.
  const G4int            fSceneId;     // Id of this instance.
  G4String               fName;
  G4int                  fViewCount;   // To determine view ids.
  G4ViewerList           fViewerList;  // Viewers.
  G4VViewer*             fpViewer;     // Current viewer.
  G4Scene*               fpScene;      // Scene for this scene handler.

  //////////////////////////////////////////////////////////////
  // Workspace...

  G4bool fReadyForTransients;  // I.e., not processing scene.
  const G4VModel*          fpModel;      // Current model.
  const G4Transform3D*     fpObjectTransformation;  // Accum'd obj. transfn.
  const G4VisAttributes*   fpVisAttribs; // Working vis attributes.
  G4int              fCurrentDepth; // Current depth of geom. hierarchy.
  G4VPhysicalVolume* fpCurrentPV;   // Current physical volume.
  G4LogicalVolume*   fpCurrentLV;   // Current logical volume.

private:

  //////////////////////////////////////////////////////////////
  // Friend function accessed only by views of this scene.

  friend void G4VViewer::ProcessView ();

  //////////////////////////////////////////////////////////////
  // Private functions, etc..

  void   ProcessScene     (G4VViewer& view);
  // Accessed by G4VViewer::ProcessView ().

};

#include "G4VSceneHandler.icc"

#endif
