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
// John Allison  19th July 1996.
//
// Class description
//
// Abstract interface class for graphics scene handlers.
// Inherits from G4VGraphicsScene, in the intercoms category, which is
// a minimal abstract interface for the GEANT4 kernel.

#ifndef G4VSCENEHANDLER_HH
#define G4VSCENEHANDLER_HH

#include "globals.hh"

#include "G4VGraphicsScene.hh"
#include "G4ViewerList.hh"
#include "G4ViewParameters.hh"
#include "G4THitsMap.hh"
#include "G4PseudoScene.hh"

class G4Scene;
class G4VGraphicsSystem;
class G4AttHolder;

class G4VSceneHandler: public G4VGraphicsScene {

  friend class G4VViewer;
  friend std::ostream& operator << (std::ostream& os, const G4VSceneHandler& s);

public: // With description

  enum MarkerSizeType {world, screen};

  G4VSceneHandler (G4VGraphicsSystem& system,
		   G4int id,
		   const G4String& name = "");

  virtual ~G4VSceneHandler ();

  ///////////////////////////////////////////////////////////////////
  // Methods for adding raw GEANT4 objects to the scene handler.  They
  // must always be called in the triplet PreAddSolid, AddSolid and
  // PostAddSolid.  The transformation and visualization attributes
  // must be set by the call to PreAddSolid.  If your graphics system
  // is sophisticated enough to handle a particular solid shape as a
  // primitive, in your derived class write a function to override one
  // or more of the following.  See the implementation of
  // G4VSceneHandler::AddSolid (const G4Box& box) for more
  // suggestions.  If not, please implement the base class invocation.

  virtual void PreAddSolid (const G4Transform3D& objectTransformation,
			    const G4VisAttributes&);
  // objectTransformation is the transformation in the world
  // coordinate system of the object about to be added, and visAttribs
  // is its visualization attributes.
  // IMPORTANT: invoke this from your polymorphic versions, e.g.:
  // void MyXXXSceneHandler::PreAddSolid
  //  (const G4Transform3D& objectTransformation,
  //   const G4VisAttributes& visAttribs) {
  //   G4VSceneHandler::PreAddSolid (objectTransformation, visAttribs);
  //   ...
  // }

  virtual void PostAddSolid ();
  // IMPORTANT: invoke this from your polymorphic versions, e.g.:
  // void MyXXXSceneHandler::PostAddSolid () {
  //   ...
  //   G4VSceneHandler::PostAddSolid ();
  // }

  // From geometry/solids/CSG
  virtual void AddSolid (const G4Box&);
  virtual void AddSolid (const G4Cons&);
  virtual void AddSolid (const G4Orb&);
  virtual void AddSolid (const G4Para&);
  virtual void AddSolid (const G4Sphere&);
  virtual void AddSolid (const G4Torus&);
  virtual void AddSolid (const G4Trap&);
  virtual void AddSolid (const G4Trd&);
  virtual void AddSolid (const G4Tubs&);

  // From geometry/solids/specific
  virtual void AddSolid (const G4Ellipsoid&);
  virtual void AddSolid (const G4Polycone&);
  virtual void AddSolid (const G4Polyhedra&);
  virtual void AddSolid (const G4TessellatedSolid&);

  // For solids not above.
  virtual void AddSolid (const G4VSolid&);

  ///////////////////////////////////////////////////////////////////
  // Methods for adding "compound" GEANT4 objects to the scene
  // handler.  These methods may either (a) invoke "user code" that
  // uses the "user interface", G4VVisManager (see, for example,
  // G4VSceneHandler, which for trajectories uses
  // G4VTrajectory::DrawTrajectory, via G4TrajectoriesModel in the
  // Modeling Category) or (b) invoke AddPrimitives below (between
  // calls to Begin/EndPrimitives) or (c) use graphics-system-specific
  // code or (d) any combination of the above.

  virtual void AddCompound (const G4VTrajectory&);
  virtual void AddCompound (const G4VHit&);
  virtual void AddCompound (const G4VDigi&);
  virtual void AddCompound (const G4THitsMap<G4double>&);
  virtual void AddCompound (const G4THitsMap<G4StatDouble>&);
  virtual void AddCompound (const G4Mesh&);

  //////////////////////////////////////////////////////////////
  // Functions for adding primitives.

  virtual void BeginModeling ();
  // IMPORTANT: invoke this from your polymorphic versions, e.g.:
  // void MyXXXSceneHandler::BeginModeling () {
  //   G4VSceneHandler::BeginModeling ();
  //   ...
  // }

  virtual void EndModeling ();
  // IMPORTANT: invoke this from your polymorphic versions, e.g.:
  // void MyXXXSceneHandler::EndModeling () {
  //   ...
  //   G4VSceneHandler::EndModeling ();
  // }

  virtual void BeginPrimitives
  (const G4Transform3D& objectTransformation = G4Transform3D());
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

  virtual void BeginPrimitives2D
  (const G4Transform3D& objectTransformation = G4Transform3D());
  // The x,y coordinates of the primitives passed to AddPrimitive are
  // intrepreted as screen coordinates, -1 < x,y < 1.  The
  // z-coordinate is ignored.
  // IMPORTANT: invoke this from your polymorphic versions, e.g.:
  // void MyXXXSceneHandler::BeginPrimitives2D
  // (const G4Transform3D& objectTransformation) {
  //   G4VSceneHandler::BeginPrimitives2D (objectTransformation);
  //   ...
  // }

  virtual void EndPrimitives2D ();
  // IMPORTANT: invoke this from your polymorphic versions, e.g.:
  // void MyXXXSceneHandler::EndPrimitives2D () {
  //   ...
  //   G4VSceneHandler::EndPrimitives2D ();
  // }

  virtual void AddPrimitive (const G4Polyline&)   = 0;
  virtual void AddPrimitive (const G4Text&)       = 0;
  virtual void AddPrimitive (const G4Circle&)     = 0;      
  virtual void AddPrimitive (const G4Square&)     = 0;      
  virtual void AddPrimitive (const G4Polymarker&);
  virtual void AddPrimitive (const G4Polyhedron&) = 0;
  virtual void AddPrimitive (const G4Plotter&);
  
  // Other virtual functions
  virtual const G4VisExtent& GetExtent() const;

  //////////////////////////////////////////////////////////////
  // Access functions.
  const G4String&     GetName           () const;
  G4int               GetSceneHandlerId () const;
  G4int               GetViewCount      () const;
  G4VGraphicsSystem*  GetGraphicsSystem () const;
  G4Scene*            GetScene          () const;
  const G4ViewerList& GetViewerList     () const;
  G4VModel*           GetModel          () const;
  G4VViewer*          GetCurrentViewer  () const;
  G4bool              GetMarkForClearingTransientStore () const;
  G4bool              IsReadyForTransients () const;
  G4bool              GetTransientsDrawnThisEvent () const;
  G4bool              GetTransientsDrawnThisRun   () const;
  const G4Transform3D& GetObjectTransformation    () const;
  void                SetName          (const G4String&);
  void          SetCurrentViewer (G4VViewer*);
  virtual void        SetScene         (G4Scene*);
  G4ViewerList& SetViewerList    ();  // Non-const so you can change.
  void          SetModel         (G4VModel*);
  void          SetMarkForClearingTransientStore (G4bool);
  // Sets flag which will cause transient store to be cleared at the
  // next call to BeginPrimitives().  Maintained by vis manager.
  void          SetTransientsDrawnThisEvent      (G4bool);
  void          SetTransientsDrawnThisRun        (G4bool);
  void          SetObjectTransformation          (const G4Transform3D&);
  // Maintained by vis manager.

  //////////////////////////////////////////////////////////////
  // Public utility functions.

  const G4Colour& GetColour ();  // To be deprecated?
  const G4Colour& GetColor  ();
  // Returns colour - checks fpVisAttribs and gets applicable colour.
  // Assumes fpVisAttribs point to the G4VisAttributes of the current object.
  // If the pointer is null, the colour is obtained from the default view
  // parameters of the current viewer.

  const G4Colour& GetColour (const G4Visible&);
  const G4Colour& GetColor  (const G4Visible&);
  // Returns colour, or viewer default colour.
  // Makes no assumptions about the validity of data member fpVisAttribs.
  // If the G4Visible has no vis attributes, i.e., the pointer is null,
  // the colour is obtained from the default view parameters of the
  // current viewer.

  const G4Colour& GetTextColour (const G4Text&);
  const G4Colour& GetTextColor  (const G4Text&);
  // Returns colour of G4Text object, or default text colour.

  G4double GetLineWidth(const G4VisAttributes*);
  // Returns line width of G4VisAttributes multiplied by GlobalLineWidthScale.

  G4ViewParameters::DrawingStyle GetDrawingStyle (const G4VisAttributes*);
  // Returns drawing style from current view parameters, unless the user
  // has forced through the vis attributes, thereby over-riding the
  // current view parameter.

  G4int GetNumberOfCloudPoints (const G4VisAttributes*) const;
  // Returns no of cloud points from current view parameters, unless the user
  // has forced through the vis attributes, thereby over-riding the
  // current view parameter.

  G4bool GetAuxEdgeVisible (const G4VisAttributes*);
  // Returns auxiliary edge visibility from current view parameters,
  // unless the user has forced through the vis attributes, thereby
  // over-riding the current view parameter.

  G4int GetNoOfSides(const G4VisAttributes*);
  // Returns no. of sides (lines segments per circle) from current
  // view parameters, unless the user has forced through the vis
  // attributes, thereby over-riding the current view parameter.

  G4double GetMarkerSize (const G4VMarker&, MarkerSizeType&);
  // Returns applicable marker size (diameter) and type (in second
  // argument).  Uses global default marker if marker sizes are not
  // set.  Multiplies by GlobalMarkerScale.

  G4double GetMarkerDiameter (const G4VMarker&, MarkerSizeType&);
  // Alias for GetMarkerSize.

  G4double GetMarkerRadius (const G4VMarker&, MarkerSizeType&);
  // GetMarkerSize / 2.

  G4ModelingParameters* CreateModelingParameters ();
  // Only the scene handler and view know what the Modeling Parameters should
  // be.  For historical reasons, the GEANT4 Visualization Environment
  // maintains its own Scene Data and View Parameters, which must be
  // converted, when needed, to Modeling Parameters.

  void DrawEvent(const G4Event*);
  // Checks scene's end-of-event model list and draws trajectories,
  // hits, etc.

  void DrawEndOfRunModels();
  // Draws end-of-run models.

  //////////////////////////////////////////////////////////////
  // Administration functions.

  template <class T> void AddSolidT (const T& solid);
  template <class T> void AddSolidWithAuxiliaryEdges (const T& solid);

  G4int IncrementViewCount ();

  virtual void ClearStore ();
  // Clears graphics database (display lists) if any.

  virtual void ClearTransientStore ();
  // Clears transient part of graphics database (display lists) if any.

  void AddViewerToList      (G4VViewer* pView);  // Add view to view List.
  void RemoveViewerFromList (G4VViewer* pView);  // Remove view from view List.

protected:

  //////////////////////////////////////////////////////////////
  // Core routine for looping over models, redrawing stored events, etc.
  // Overload with care (see, for example,
  // G4OpenGLScenehandler::ProcessScene).
  virtual void ProcessScene ();

  //////////////////////////////////////////////////////////////
  // Default routine used by default AddSolid ().
  virtual void RequestPrimitives (const G4VSolid& solid);

  //////////////////////////////////////////////////////////////
  // Other internal routines...

  virtual G4DisplacedSolid* CreateSectionSolid ();
  virtual G4DisplacedSolid* CreateCutawaySolid ();
  // Generic clipping using the BooleanProcessor in graphics_reps is
  // implemented in this class.  Subclasses that implement their own
  // clipping should provide an override that returns zero.

  void LoadAtts(const G4Visible&, G4AttHolder*);
  // Load G4AttValues and G4AttDefs associated with the G4Visible
  // object onto the G4AttHolder object.  It checks fpModel, and also
  // loads the G4AttValues and G4AttDefs from G4PhysicalVolumeModel,
  // G4VTrajectory, G4VTrajectoryPoint, G4VHit or G4VDigi, as
  // appropriate.  The G4AttHolder object is an object of a class that
  // publicly inherits G4AttHolder - see, e.g., SoG4Polyhedron in the
  // Open Inventor driver.  G4AttHolder deletes G4AttValues in its
  // destructor to ensure proper clean-up of G4AttValues.

  //////////////////////////////////////////////////////////////
  // Special mesh rendering utilities...

  struct NameAndVisAtts {
    NameAndVisAtts(const G4String& name = "",const G4VisAttributes& visAtts = G4VisAttributes())
    : fName(name),fVisAtts(visAtts) {}
    G4String fName;
    G4VisAttributes fVisAtts;
  };

  class PseudoSceneFor3DRectMeshPositions: public G4PseudoScene {
  public:
    PseudoSceneFor3DRectMeshPositions
    (G4PhysicalVolumeModel* pvModel  // input
     , G4int depth // input...the following are outputs by reference
     , std::multimap<const G4Material*,const G4ThreeVector>& positionByMaterial
     , std::map<const G4Material*,NameAndVisAtts>& nameAndVisAttsByMaterial)
    : fpPVModel(pvModel)
    , fDepth(depth)
    , fPositionByMaterial(positionByMaterial)
    , fNameAndVisAttsByMaterial(nameAndVisAttsByMaterial)
    {}
  private:
    using G4PseudoScene::AddSolid;  // except for...
    void AddSolid(const G4Box&) override;
    void ProcessVolume(const G4VSolid&) override {
      // Do nothing if uninteresting solids found, e.g., the container if not marked invisible.
    }
    G4PhysicalVolumeModel* fpPVModel;
    G4int fDepth;
    std::multimap<const G4Material*,const G4ThreeVector>& fPositionByMaterial;
    std::map<const G4Material*,NameAndVisAtts>& fNameAndVisAttsByMaterial;
  };

  class PseudoSceneForTetVertices: public G4PseudoScene {
  public:
    PseudoSceneForTetVertices
    (G4PhysicalVolumeModel* pvModel  // input
     , G4int depth // input...the following are outputs by reference
     , std::multimap<const G4Material*,std::vector<G4ThreeVector>>& verticesByMaterial
     , std::map<const G4Material*,NameAndVisAtts>& nameAndVisAttsByMaterial)
    : fpPVModel(pvModel)
    , fDepth(depth)
    , fVerticesByMaterial(verticesByMaterial)
    , fNameAndVisAttsByMaterial(nameAndVisAttsByMaterial)
    {}
  private:
    using G4PseudoScene::AddSolid;  // except for...
    void AddSolid(const G4VSolid& solid) override;
    void ProcessVolume(const G4VSolid&) override {
      // Do nothing if uninteresting solids found, e.g., the container if not marked invisible.
    }
    G4PhysicalVolumeModel* fpPVModel;
    G4int fDepth;
    std::multimap<const G4Material*,std::vector<G4ThreeVector>>& fVerticesByMaterial;
    std::map<const G4Material*,NameAndVisAtts>& fNameAndVisAttsByMaterial;
  };

  void StandardSpecialMeshRendering(const G4Mesh&);
  // Standard way of special mesh rendering.
  // MySceneHandler::AddCompound(const G4Mesh& mesh) may use this if
  // appropriate or implement its own special mesh rendereing.

  void Draw3DRectMeshAsDots(const G4Mesh&);
  // For a rectangular 3-D mesh, draw as coloured dots by colour and material,
  // one dot randomly placed in each visible mesh cell.

  void Draw3DRectMeshAsSurfaces(const G4Mesh&);
  // For a rectangular 3-D mesh, draw as surfaces by colour and material
  // with inner shared faces removed.

  void DrawTetMeshAsDots(const G4Mesh&);
  // For a tetrahedron mesh, draw as coloured dots by colour and material,
  // one dot randomly placed in each visible mesh cell.

  void DrawTetMeshAsSurfaces(const G4Mesh&);
  // For a tetrahedron mesh, draw as surfaces by colour and material
  // with inner shared faces removed.

  G4ThreeVector GetPointInBox(const G4ThreeVector& pos,
                              G4double halfX,
                              G4double halfY,
                              G4double halfZ) const;
  // Sample a random point inside the box

  G4ThreeVector GetPointInTet(const std::vector<G4ThreeVector>& vertices) const;
  // Sample a random point inside the tetrahedron

  //////////////////////////////////////////////////////////////
  // Data members

  G4VGraphicsSystem& fSystem;          // Graphics system.
  const G4int        fSceneHandlerId;  // Id of this instance.
  G4String           fName;
  G4int              fViewCount;       // To determine view ids.
  G4ViewerList       fViewerList;      // Viewers.
  G4VViewer*         fpViewer;         // Current viewer.
  G4Scene*           fpScene;          // Scene for this scene handler.
  G4bool             fMarkForClearingTransientStore;
  G4bool             fReadyForTransients;  // I.e., not processing the
			                   // run-duration part of scene.
  G4bool             fTransientsDrawnThisEvent;  // Maintained by vis
  G4bool             fTransientsDrawnThisRun;    // manager.
  G4bool             fProcessingSolid; // True if within Pre/PostAddSolid.
  G4bool             fProcessing2D;    // True for 2D.
  G4VModel*          fpModel;          // Current model.
  G4Transform3D fObjectTransformation; // Current accumulated
				       // object transformation.
  G4int              fNestingDepth;    // For Begin/EndPrimitives.
  const G4VisAttributes* fpVisAttribs; // Working vis attributes.
  const G4Transform3D fIdentityTransformation;

private:

  G4VSceneHandler (const G4VSceneHandler&);
  G4VSceneHandler& operator = (const G4VSceneHandler&);
};

#include "G4VSceneHandler.icc"

#endif
