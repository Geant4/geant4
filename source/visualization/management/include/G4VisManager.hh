// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisManager.hh,v 1.5 1999-02-07 17:32:15 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

// The GEANT4 Visualization Manager - John Allison 02/Jan/1996.

// This is a "Singleton", i.e., only one instance may exist.  A
// G4Exception is thrown if an attempt is made to instantiate more
// than one.

// Graphics system registration is normally done through the protected
// pure virtual function RegisterGraphicsSystems called from
// Initialise ().  You must define your own subclass and implement
// RegisterGraphicsSystems - for an example see
// visualization/include/MyVisManager.hh/cc.  You can also use the
// public function RegisterGraphicsSystem (new MyGraphicsSystem) if
// you have your own graphics system.

// The VisManager creates graphics systems, scenes, scene handlers and
// viewers and manages them.  You can have any number.  It has the
// concept of a "current viewer", and the "current scene handler", the
// "current scene" and the "current graphics system" which go with it.
// You can select the current viewer.  Most of the the operations of
// the VisManager take place with the current viewer, in particular,
// the Draw operations.

// Each scene comprises drawable objects such as detector components
// and hits when appropriate.  A scene handler translates a scene into
// graphics-system-specific function calls and, possibly, a
// graphics-system-dependent database - display lists, scene graphs,
// etc.  Each viewer has its "view parameters".

// G4VisManager is "state dependent", i.e., it is notified on change
// of state (G4ApplicationState).  This is used to draw hits and
// trajectories in the current scene at the end of event, as required.

#ifndef G4VISMANAGER_HH
#define G4VISMANAGER_HH

#include "G4VVisManager.hh"

#include "globals.hh"
#include "G4GraphicsSystemList.hh"
#include "G4SceneHandlerList.hh"
#include "G4Scene.hh"
#include "G4SceneList.hh"
#include "G4ViewParameters.hh"
#include "G4Transform3D.hh"
#include "G4UImessenger.hh"
#include "G4VStateDependent.hh"

#include <rw/tpordvec.h>

class G4VisManMessenger;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSolid;
class G4VGraphicsSystem;
class G4VSceneHandler;
class G4VViewer;
class G4Polyline;
class G4Text;
class G4Circle;
class G4Square;
class G4Polymarker;
class G4Polyhedron;
class G4NURBS;

class G4VisManager: public G4VVisManager, public G4VStateDependent {

  // List of classes and functions which need access to private
  // members of G4VisManager.  This is mainly to obtain access to
  // GetInstance (), which is private.  The correct way for normal
  // users to obtain a pointer to the vis manager is with
  // G4VVisManager::GetConcreteInstance (), always testing for
  // non-zero.

  friend void G4OpenGLXmViewerSecondaryLoopPostAction ();  // Mmmm!

  friend class G4VSceneHandler;
  friend class G4VViewer;

  friend class ostream & operator <<
  (class ostream &, const class G4VGraphicsSystem &);

  friend class ostream & operator <<
  (class ostream &, const class G4VSceneHandler &);

  // Now classes associated with the old commands...
  friend class G4VisManMessenger;
  friend class G4VisCommandCameraReset;
  friend class G4VisCommandClearScene;
  friend class G4VisCommandClearView;
  friend class G4VisCommandClearViewAndScene;
  friend class G4VisCommandCopyAll;
  friend class G4VisCommandCopyScene;
  friend class G4VisCommandCopyView;
  friend class G4VisCommandCreateViewNewScene;
  friend class G4VisCommandCreateViewNewView;
  friend class G4VisCommandDeleteScene;
  friend class G4VisCommandDeleteView;
  friend class G4VisCommandDrawCurrent;
  friend class G4VisCommandLightsMoveWithCamera;
  friend class G4VisCommandRefreshView;
  friend class G4VisCommandSetCulling;
  friend class G4VisCommandSetCullCoveredDaughters;
  friend class G4VisCommandSetCullInvisible;
  friend class G4VisCommandShowView;

public:

  G4VisManager ();

  virtual ~G4VisManager ();

private:

  static G4VisManager* GetInstance ();
  // Returns pointer to itself.  Throws a G4Exception if called before
  // instantiation.  Private so that only friends can use; the normal
  // user should instead use G4VVisManager::GetConcreteInstance () to
  // get a "higher level" pointer for general use - but always test
  // for non-zero.

public:

  void Initialise ();
  void Initialize ();  // Alias Initialise ().
  G4bool RegisterGraphicsSystem (G4VGraphicsSystem* pSystem);

  //////////////////////////////////////////////////////////////////////
  // Drawing routines!

  void ClearView ();
  // Clear visible window of current viewer (both buffers of a double
  // buffered system).

  void Draw ();
  // Draw current scene in current view.

  void Show ();
  // Show current view (initiates post-procesing of graphical
  // databases for graphics systems which require it.)

  // ///////////////////////////////////////////////////////////////
  // Now functions for drawing various visualization primitives,
  // useful for representing hits, digis, etc.

  void Draw (const G4Polyline&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity);

  void Draw (const G4Text&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity);

  void Draw (const G4Circle&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity);

  void Draw (const G4Square&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity);

  void Draw (const G4Polymarker&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity);

  void Draw (const G4Polyhedron&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity);

  void Draw (const G4NURBS&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity);

  // //////////////////////////////////////////////////////////////////
  // Now functions for drawing a GEANT4 geometry object.  Note that
  // the 2nd argument overrides any visualization attributes that are
  // associated with the object itself.

  void Draw (const G4VSolid&, const G4VisAttributes&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity);

  void Draw (const G4LogicalVolume&, const G4VisAttributes&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity);

  void Draw (const G4VPhysicalVolume&, const G4VisAttributes&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity);

  G4bool Notify (G4ApplicationState requestedState);
  // This is called on change of state (G4ApplicationState).  It is
  // used to draw hits and trajectories if included in the current
  // scene at the end of event, as required.

  ////////////////////////////////////////////////////////////////////////
  // Administration routines.

  void CopyViewParameters ();
  // Copy view parameters of current viewer into current view parameters.

  void CreateSceneHandler (G4String name = "");
  // Creates scene handler for the current system.

  void CreateViewer  (G4String name = "");
  // Creates viewer for the current scene handler.

  void DeleteCurrentSceneHandler ();
  // Leaves current scene handler and viewer undefined!

  void DeleteCurrentViewer ();
  // Leaves current scene and view undefined!

  void GeometryHasChanged ();
  // Used by run manager to notify change.

  void RefreshCurrentView  ();
  // Soft clear, then redraw.

  // These can go when OLD STYLE commands go...
  void PrintCurrentSystem  () const;
  void PrintCurrentSystems () const;
  void PrintCurrentScene   () const;
  void PrintCurrentView    () const;

  /////////////////////////////////////////////////////////////////////
  // Access functions.

  G4VGraphicsSystem*           GetCurrentGraphicsSystem    () const;
  G4Scene*                     GetCurrentScene             () const;
  G4VSceneHandler*             GetCurrentSceneHandler      () const;
  G4VViewer*                   GetCurrentViewer            () const;
  const G4ViewParameters&      GetCurrentViewParameters    () const;
  const G4GraphicsSystemList&  GetAvailableGraphicsSystems ();
  // The above is non-const because it checks and updates the List by
  // calling RegisterGraphicsSystems() if no graphics systems are
  // already registered.
  const G4SceneHandlerList&    GetAvailableSceneHandlers   () const;
  const G4SceneList&           GetSceneList                () const;
  G4int                        GetVerboseLevel             () const;
  G4bool                       IsValidView                 ();
  void              SetCurrentGraphicsSystemAndCreateViewer
                                                (G4VGraphicsSystem* pSystem);
  void              SetCurrentGraphicsSystem    (G4VGraphicsSystem* pSystem);
  void              SetCurrentScene             (G4Scene*);
  void              SetCurrentSceneHandler      (G4VSceneHandler* pScene);
  void              SetCurrentViewer            (G4VViewer* pView);
  G4ViewParameters& SetCurrentViewParameters    ();  // Returns lvalue.
  G4SceneHandlerList& SetAvailableSceneHandlers ();  // Returns lvalue.
  G4SceneList&      SetSceneList                ();  // Returns lvalue.
  void              SetVerboseLevel             (G4int vLevel);

protected:

  virtual void RegisterGraphicsSystems () = 0;
  void RegisterMessengers              ();   // Command messengers.
  void PrintAllGraphicsSystems         () const;
  void PrintInstalledGraphicsSystems   () const;
  void PrintAvailableGraphicsSystems   () const;
  void PrintInvalidPointers            () const;
  static G4VisManager*  fpInstance;         // Pointer to single instance.
  G4bool                fInitialised;
  G4VGraphicsSystem*    fpGraphicsSystem;   // Current graphics system.
  G4Scene*              fpScene;            // Current scene.
  G4VSceneHandler*      fpSceneHandler;     // Current scene handler.
  G4VViewer*            fpViewer;           // Current viewer.
  G4GraphicsSystemList  fAvailableGraphicsSystems;
  G4SceneList           fSceneList;
  G4SceneHandlerList    fAvailableSceneHandlers;
  G4ViewParameters      fVP;                // Current viewing parameters.
  G4int                 fVerbose;           // Verbosity level 0-10.
  G4VisManMessenger*    fpMessenger;        // Pointer to messenger.
  RWTPtrOrderedVector <G4UImessenger> fMessengerList;
};

#include "G4VisManager.icc"

#endif
