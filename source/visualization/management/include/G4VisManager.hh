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
// $Id: G4VisManager.hh 102802 2017-02-22 15:19:00Z gcosmo $
//
// 

// Class Description:
//
// The GEANT4 Visualization Manager - John Allison 02/Jan/1996.
//
// G4VisManager is a "Singleton", i.e., only one instance of it or any
// derived class may exist.  A G4Exception is thrown if an attempt is
// made to instantiate more than one.
//
// It is also an abstract class, so the user must derive his/her own
// class from G4VisManager, implement the pure virtual function
// RegisterGraphicsSystems, and instantiate an object of the derived
// class - for an example see
// visualization/include/G4VisExecutive.hh/icc.
//
// The recommended way for users to obtain a pointer to the vis
// manager is with G4VVisManager::GetConcreteInstance (), being always
// careful to test for non-zero.  This pointer is non-zero only when
// (a) an object of the derived class exists and (b) when there is a
// valid viewer available.
//
// Graphics system registration is normally done through the protected
// pure virtual function RegisterGraphicsSystems called from
// Initialise ().  You can also use the public function
// RegisterGraphicsSystem (new MyGraphicsSystem) if you have your own
// graphics system. A graphics system is, in effect, a factory for
// scene handlers and viewers.
//
// The VisManager creates and manages graphics systems, scenes, scene
// handlers, viewers and some models and model makers.  You can have
// any number.  It has the concept of a "current viewer", and the
// "current scene handler", the "current scene" and the "current
// graphics system" which go with it.  You can select the current
// viewer.  Most of the the operations of the VisManager take place
// with the current viewer, in particular, the Draw operations.
//
// Each scene comprises drawable objects such as detector components
// and trajectories, hits and digis when appropriate.  A scene handler
// translates a scene into graphics-system-specific function calls
// and, possibly, a graphics-system-dependent database - display
// lists, scene graphs, etc.  Each viewer has its "view parameters"
// (see class description of G4ViewParameters for available parameters
// and also for a description of the concept of a "standard view" and
// all that).
//
// A friend class G4VisStateDependent is "state dependent", i.e., it
// is notified on change of state (G4ApplicationState).  This is used
// to message the G4VisManager to draw hits, digis and trajectories in
// the current scene at the end of event, as required.

#ifndef G4VISMANAGER_HH
#define G4VISMANAGER_HH

// Temporary definition until Xeon Phi can handle full C++11.
#ifndef __MIC__
#define G4VIS_USE_STD11
#endif

#include "G4VVisManager.hh"

#include "globals.hh"
#include "G4GraphicsSystemList.hh"
#include "G4ModelingParameters.hh"
#include "G4NullModel.hh"
#include "G4SceneHandlerList.hh"
#include "G4SceneList.hh"
#include "G4TrajectoriesModel.hh"
#include "G4Transform3D.hh"
#include "G4UImessenger.hh"

#include <iostream>
#include <vector>
#include <map>

#include "G4Threading.hh"

class G4Scene;
class G4UIcommand;
class G4UImessenger;
class G4VisStateDependent;
class G4VTrajectoryModel;
class G4VUserVisAction;
template <typename> class G4VFilter;
template <typename> class G4VisFilterManager;
template <typename> class G4VisModelManager;
template <typename> class G4VModelFactory;
class G4Event;

// Useful typedef's
typedef G4VModelFactory<G4VTrajectoryModel> G4TrajDrawModelFactory;
typedef G4VModelFactory<G4VFilter<G4VTrajectory> > G4TrajFilterFactory;
typedef G4VModelFactory<G4VFilter<G4VHit> > G4HitFilterFactory;
typedef G4VModelFactory<G4VFilter<G4VDigi> > G4DigiFilterFactory;

class G4VisManager: public G4VVisManager {

  // Management friends...
  friend class G4VSceneHandler;
  friend class G4VViewer;
  friend class G4VisStateDependent;
  friend class G4VisCommandList;

  // operator << friends...
  friend std::ostream& operator << (std::ostream&, const G4VGraphicsSystem&);
  friend std::ostream& operator << (std::ostream&, const G4VSceneHandler&);

public: // With description

  enum Verbosity {
    quiet,         // Nothing is printed.
    startup,       // Startup and endup messages are printed...
    errors,        // ...and errors...
    warnings,      // ...and warnings...
    confirmations, // ...and confirming messages...
    parameters,    // ...and parameters of scenes and views...
    all            // ...and everything available.
  };
  // Simple graded message scheme.

protected: // With description

  G4VisManager (const G4String& verbosityString = "warnings");
  // The constructor is protected so that an object of the derived
  // class may be constructed.

public: // With description

  virtual ~G4VisManager ();

private:

  // Private copy constructor and assigment operator - copying and
  // assignment not allowed.  Keeps CodeWizard happy.
  G4VisManager (const G4VisManager&);
  G4VisManager& operator = (const G4VisManager&);

public:
  static G4VisManager* GetInstance ();
  // Returns pointer to itself.  Throws a G4Exception if called before
  // instantiation.  Intended only for use within the vis category; the
  // normal user should instead use G4VVisManager::GetConcreteInstance()
  // to get a "higher level" pointer for general use - but always test
  // for non-zero.

public: // With description

  void Initialise ();
  void Initialize ();  // Alias Initialise ().

  // Optional registration of user vis actions.  Added to scene with
  // /vis/scene/add/userAction.
  void RegisterRunDurationUserVisAction
  (const G4String& name, G4VUserVisAction*,
   const G4VisExtent& = G4VisExtent());
  void RegisterEndOfEventUserVisAction
  (const G4String& name, G4VUserVisAction*,
   const G4VisExtent& = G4VisExtent());
  void RegisterEndOfRunUserVisAction
  (const G4String& name, G4VUserVisAction*,
   const G4VisExtent& = G4VisExtent());
  void SetUserAction
  (G4VUserVisAction* pVisAction,
   const G4VisExtent& = G4VisExtent());
  // SetUserAction is deprecated.  Use RegisterRunDurationUserVisAction
  // or other of the above.
  void SetUserActionExtent (const G4VisExtent&);  //Legacy: deprecated.


  G4bool RegisterGraphicsSystem (G4VGraphicsSystem*);
  // Register an individual graphics system.  Normally this is done in
  // a sub-class implementation of the protected virtual function,
  // RegisterGraphicsSystems.  See, e.g., G4VisExecutive.icc.

  void RegisterModelFactory(G4TrajDrawModelFactory* factory);
  // Register trajectory draw model factory. Assumes ownership of factory.

  void RegisterModel(G4VTrajectoryModel* model);
  // Register trajectory model. Assumes ownership of model.

  void RegisterModelFactory(G4TrajFilterFactory* factory);
  // Register trajectory filter model factory. Assumes ownership of factory.

  void RegisterModel(G4VFilter<G4VTrajectory>* filter);
  // Register trajectory filter model. Assumes ownership of model.

  void RegisterModelFactory(G4HitFilterFactory* factory);
  // Register trajectory hit model factory. Assumes ownership of factory.

  void RegisterModel(G4VFilter<G4VHit>* filter);
  // Register trajectory hit model. Assumes ownership of model.

  void RegisterModelFactory(G4DigiFilterFactory* factory);
  // Register trajectory digi model factory. Assumes ownership of factory.

  void RegisterModel(G4VFilter<G4VDigi>* filter);
  // Register trajectory digi model. Assumes ownership of model.

  void SelectTrajectoryModel(const G4String& model);
  // Set default trajectory model. Useful for use in compiled code

  void RegisterMessenger(G4UImessenger* messenger);
  // Register messenger. Assumes ownership of messenger.

  /////////////////////////////////////////////////////////////////
  // Now functions that implement the pure virtual functions of
  // G4VVisManager for drawing various visualization primitives, useful
  // for representing hits, digis, etc.

  void Draw (const G4Circle&,
    const G4Transform3D& objectTransformation = G4Transform3D());

  void Draw (const G4Polyhedron&,
    const G4Transform3D& objectTransformation = G4Transform3D());

  void Draw (const G4Polyline&,
    const G4Transform3D& objectTransformation = G4Transform3D());

  void Draw (const G4Polymarker&,
    const G4Transform3D& objectTransformation = G4Transform3D());

  void Draw (const G4Scale&,
    const G4Transform3D& objectTransformation = G4Transform3D());

  void Draw (const G4Square&,
    const G4Transform3D& objectTransformation = G4Transform3D());

  void Draw (const G4Text&,
    const G4Transform3D& objectTransformation = G4Transform3D());

  void Draw2D (const G4Circle&,
    const G4Transform3D& objectTransformation = G4Transform3D());

  void Draw2D (const G4Polyhedron&,
    const G4Transform3D& objectTransformation = G4Transform3D());

  void Draw2D (const G4Polyline&,
    const G4Transform3D& objectTransformation = G4Transform3D());

  void Draw2D (const G4Polymarker&,
    const G4Transform3D& objectTransformation = G4Transform3D());

  void Draw2D (const G4Square&,
    const G4Transform3D& objectTransformation = G4Transform3D());

  void Draw2D (const G4Text&,
    const G4Transform3D& objectTransformation = G4Transform3D());

  ////////////////////////////////////////////////////////////////////
  // Now functions that implement the pure virtual functions of
  // G4VVisManager for drawing a GEANT4 object.  Note that the
  // visualization attributes needed in some cases override any
  // visualization attributes that are associated with the object
  // itself - thus you can, for example, change the colour of a
  // physical volume.

  void Draw (const G4VTrajectory&);

  void Draw (const G4VHit&);

  void Draw (const G4VDigi&);

  void Draw (const G4LogicalVolume&, const G4VisAttributes&,
    const G4Transform3D& objectTransformation = G4Transform3D());

  void Draw (const G4VPhysicalVolume&, const G4VisAttributes&,
    const G4Transform3D& objectTransformation = G4Transform3D());

  void Draw (const G4VSolid&, const G4VisAttributes&,
    const G4Transform3D& objectTransformation = G4Transform3D());

  //////////////////////////////////////////////////////////////////////
  // Optional methods that you may use to bracket a series of Draw
  // messages that have identical objectTransformation to improve
  // drawing speed.  Use Begin/EndDraw for a series of Draw messages,
  // Begin/EndDraw2D for a series of Draw2D messages.  Do not mix Draw
  // and Draw2D messages.

  void BeginDraw
  (const G4Transform3D& objectTransformation = G4Transform3D());

  void EndDraw ();

  void BeginDraw2D
  (const G4Transform3D& objectTransformation = G4Transform3D());

  void EndDraw2D ();

  ////////////////////////////////////////////////////////////////////////
  // Now other pure virtual functions of G4VVisManager...

  void GeometryHasChanged ();
  // Used by run manager to notify change.

  void IgnoreStateChanges(G4bool);
  // This method shoud be invoked by a class that has its own event loop,
  // such as the RayTracer, material scanner, etc. If the argument is true,
  // the following state changes among Idle, GeomClosed and EventProc are
  // caused by such a class, and thus not by the ordinary event simulation.
  // The same method with false should be invoked once such an event loop
  // is over.

  void NotifyHandlers();
  // Notify scene handlers (G4VGraphicsScene objects) that the scene
  // has changed so that they may rebuild their graphics database, if
  // any, and redraw all views.

  void DispatchToModel(const G4VTrajectory&);
  // Draw the trajectory.

  G4bool FilterTrajectory(const G4VTrajectory&);
  G4bool FilterHit(const G4VHit&);
  G4bool FilterDigi(const G4VDigi&);

#ifdef G4MULTITHREADED

  virtual void SetUpForAThread();
  // This method is invoked by G4WorkerRunManager

  static G4ThreadFunReturnType G4VisSubThread(G4ThreadFunArgType);
  // Vis sub-thread function.

#endif

  ////////////////////////////////////////////////////////////////////////
  // Administration routines.

  void CreateSceneHandler (const G4String& name = "");
  // Creates scene handler for the current system.

  void CreateViewer (const G4String& name = "", const G4String& XGeometry = "");
  // Creates viewer for the current scene handler.

private:

  void BeginOfRun ();

  void BeginOfEvent ();

  void EndOfEvent ();
  // This is called on change of state (G4ApplicationState).  It is
  // used to draw hits, digis and trajectories if included in the
  // current scene at the end of event, as required.

  void EndOfRun ();

public: // With description

  /////////////////////////////////////////////////////////////////////
  // Access functions.

  void Enable();
  void Disable();
  // Global enable/disable functions.

  const G4VTrajectoryModel* CurrentTrajDrawModel() const;

  struct UserVisAction {
    UserVisAction(const G4String& name, G4VUserVisAction* pUserVisAction)
      :fName(name), fpUserVisAction(pUserVisAction) {}
    G4String fName;
    G4VUserVisAction* fpUserVisAction;
  };
  const std::vector<UserVisAction>& GetRunDurationUserVisActions () const;
  const std::vector<UserVisAction>& GetEndOfEventUserVisActions  () const;
  const std::vector<UserVisAction>& GetEndOfRunUserVisActions    () const;
  const std::map<G4VUserVisAction*,G4VisExtent>& GetUserVisActionExtents () const;
  G4VGraphicsSystem*           GetCurrentGraphicsSystem    () const;
  G4Scene*                     GetCurrentScene             () const;
  G4VSceneHandler*             GetCurrentSceneHandler      () const;
  G4VViewer*                   GetCurrentViewer            () const;
  const G4GraphicsSystemList&  GetAvailableGraphicsSystems ();
  // The above is non-const because it checks and updates the List by
  // calling RegisterGraphicsSystems() if no graphics systems are
  // already registered.
  const G4SceneHandlerList&    GetAvailableSceneHandlers   () const;
  const G4SceneList&           GetSceneList                () const;
  static Verbosity             GetVerbosity                ();
  G4bool                       GetTransientsDrawnThisRun   () const;
  G4bool                       GetTransientsDrawnThisEvent () const;
  const G4Event*               GetRequestedEvent           () const;
  G4bool                       GetAbortReviewKeptEvents    () const;
  const G4ViewParameters&      GetDefaultViewParameters    () const;
#ifdef G4MULTITHREADED
  G4int                        GetMaxEventQueueSize        () const;
  G4bool                       GetWaitOnEventQueueFull     () const;
#endif

  void              SetCurrentGraphicsSystem    (G4VGraphicsSystem*);
  void              SetCurrentScene             (G4Scene*);
  void              SetCurrentSceneHandler      (G4VSceneHandler*);
  void              SetCurrentViewer            (G4VViewer*);
  G4SceneHandlerList& SetAvailableSceneHandlers ();  // Returns lvalue.
  G4SceneList&      SetSceneList                ();  // Returns lvalue.
  void              SetVerboseLevel             (G4int);
  void              SetVerboseLevel             (const G4String&);
  void              SetVerboseLevel             (Verbosity);
  void              SetEventRefreshing          (G4bool);
  void              ResetTransientsDrawnFlags   ();
  void              SetTransientsDrawnThisRun   (G4bool);
  void              SetTransientsDrawnThisEvent (G4bool);
  // If non-zero, requested event is used in G4VSceneHandler::ProcessScene.
  void              SetRequestedEvent           (const G4Event*);
  void              SetAbortReviewKeptEvents    (G4bool);
  void              SetDefaultViewParameters    (const G4ViewParameters&);
#ifdef G4MULTITHREADED
  void              SetMaxEventQueueSize        (G4int);
  void              SetWaitOnEventQueueFull     (G4bool);
#endif

  /////////////////////////////////////////////////////////////////////
  // Utility functions.

  G4String ViewerShortName (const G4String& viewerName) const;
  // Returns shortened version of viewer name, i.e., up to first space,
  // if any.

  G4VViewer* GetViewer (const G4String& viewerName) const;
  // Returns zero if not found.  Can use long or short name, but find
  // is done on short name.

  static Verbosity GetVerbosityValue(const G4String&);
  // Returns verbosity given a string.  (Uses first character only.)

  static Verbosity GetVerbosityValue(G4int);
  // Returns verbosity given an integer.  If integer is out of range,
  // selects verbosity at extreme of range.

  static G4String VerbosityString(Verbosity);
  // Converts the verbosity into a string for suitable for printing.
  
  static std::vector<G4String> VerbosityGuidanceStrings;
  // Guidance on the use of visualization verbosity.

protected:

  virtual void RegisterGraphicsSystems () = 0;
  // The sub-class must implement and make successive calls to
  // RegisterGraphicsSystem.

  virtual void RegisterModelFactories();
  // Sub-class must register desired models

  void RegisterMessengers              ();   // Command messengers.

  const G4int           fVerbose;
  // fVerbose is kept for backwards compatibility for some user
  // examples.  (It is used in the derived user vis managers to print
  // available graphics systems.)  It is initialised to 1 in the
  // constructor and cannot be changed.

  void PrintAvailableGraphicsSystems   (Verbosity) const;

private:

  // Function templates to implement the Draw methods (to avoid source
  // code duplication).
  template <class T> void DrawT
  (const T& graphics_primitive, const G4Transform3D& objectTransform);
  template <class T> void DrawT2D
  (const T& graphics_primitive, const G4Transform3D& objectTransform);

  void PrintAvailableModels            (Verbosity) const;
  void PrintAvailableColours           (Verbosity) const;
  void PrintAvailableUserVisActions   (Verbosity) const;
  void PrintInvalidPointers            () const;
  G4bool IsValidView ();
  // True if view is valid.  Prints messages and sanitises various data.
  void ClearTransientStoreIfMarked();
  // Clears transient store of current scene handler if it is marked
  // for clearing.  Assumes view is valid.

  static G4VisManager*  fpInstance;         // Pointer to single instance. 
  G4bool                fInitialised;
  std::vector<UserVisAction> fRunDurationUserVisActions;
  std::vector<UserVisAction> fEndOfEventUserVisActions;
  std::vector<UserVisAction> fEndOfRunUserVisActions;
  std::map<G4VUserVisAction*,G4VisExtent> fUserVisActionExtents;
  G4VGraphicsSystem*    fpGraphicsSystem;   // Current graphics system.
  G4Scene*              fpScene;            // Current scene.
  G4VSceneHandler*      fpSceneHandler;     // Current scene handler.
  G4VViewer*            fpViewer;           // Current viewer.
  G4GraphicsSystemList  fAvailableGraphicsSystems;
  G4SceneList           fSceneList;
  G4SceneHandlerList    fAvailableSceneHandlers;
  static Verbosity            fVerbosity;
  std::vector<G4UImessenger*> fMessengerList;
  std::vector<G4UIcommand*>   fDirectoryList;
  G4VisStateDependent*  fpStateDependent;   // Friend state dependent class.
  G4bool                fEventRefreshing;
  G4bool                fTransientsDrawnThisRun;
  G4bool                fTransientsDrawnThisEvent;
  G4int                 fNoOfEventsDrawnThisRun;
  G4int                 fNKeepRequests;
  G4bool                fEventKeepingSuspended;
  G4bool                fKeptLastEvent;
  const G4Event*        fpRequestedEvent; // If non-zero, scene handler uses.
  G4bool                fAbortReviewKeptEvents;
  G4ViewParameters      fDefaultViewParameters;
  G4bool                fIsDrawGroup;
  G4int                 fDrawGroupNestingDepth;
  G4bool                fIgnoreStateChanges;
#ifdef G4MULTITHREADED
  G4int                 fMaxEventQueueSize;
  G4bool                fWaitOnEventQueueFull;
#endif

  // Trajectory draw model manager
  G4VisModelManager<G4VTrajectoryModel>* fpTrajDrawModelMgr;
  
  // Trajectory filter model manager
  G4VisFilterManager<G4VTrajectory>* fpTrajFilterMgr;

  // Hit filter model manager
  G4VisFilterManager<G4VHit>* fpHitFilterMgr;

  // Digi filter model manager
  G4VisFilterManager<G4VDigi>* fpDigiFilterMgr;
};

#include "G4VisManager.icc"

#endif
