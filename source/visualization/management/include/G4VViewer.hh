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
// John Allison  27th March 1996
//
// Class description
//
// Abstract interface class for graphics viewers.

#ifndef G4VVIEWER_HH
#define G4VVIEWER_HH

#include "globals.hh"

#include "G4ViewParameters.hh"
#include "G4PhysicalVolumeModel.hh"

class G4VSceneHandler;

class G4VViewer {

public: // With description

  friend std::ostream& operator << (std::ostream& os, const G4VViewer& v);

  G4VViewer (G4VSceneHandler&, G4int id, const G4String& name = "");
  virtual ~G4VViewer ();

  virtual void Initialise ();
  // Called immediately after construction for those operations that
  // must await complete contruction of viewer and all its bases.  For
  // example, if this class (G4VViewer) is inherited virtually, as in
  // the OpenGL sub-category, it will not be fully constructed until
  // *after* the the derived viewer (this is the rule about order of
  // construction for virtual inheritance), so the derived viewer may
  // not use information in G4VViewer in its contructor.  Hence such
  // code must be in Initialise().

  //////////////////////////////////////////////////////////////
  // View manipulation functions.

  virtual void ResetView ();
  // Reset view parameters to default, including sub-class parameters, if any.
  // The sub-class should always invoke the base class implementation, i.e:
  // virtual void SubClass::ResetView () {
  //   G4VViewer::ResetView();
  //   // Then reset sub-class parameters
  //   ...

  virtual void SetView () = 0;
  // Take view parameters and work out model/view transformation,
  // projection transformation, lighting, etc.

  virtual void ClearView () = 0;
  // Clear screen/viewing buffers.

  virtual void DrawView () = 0;
  // Draw view of the scene currently attached to the scene handler -
  // see example of a minimal function at end of this file.
  
  void RefreshView ();
  // Simply invokes SetView, ClearView, DrawView.

  virtual void ShowView ();
  // Show view (for graphics systems which require to process
  // all drawn objects before finalising the view).

  virtual void FinishView ();
  // Called at the end of drawing scene.  Used to flush streams, or
  // swap buffers.  (Perhaps it is inappropriately named, perhaps its
  // function could be incorporated into EndModeling ().  It marks the
  // end of scene drawing; be aware hits and digi drawing may Follow.
  // It is not yet the end of all drawing; that is signalled by
  // ShowView ().)

  std::vector<G4ThreeVector> ComputeFlyThrough(G4Vector3D*);

#ifdef G4MULTITHREADED
  // Note: the order of calling of MovingToVisSubThread and SwitchToVisSubThread
  // is undefined, so you may need to implement mutexes to ensure your preferred
  // order - see, e.g., G4OpenGLQtViewer. To summarise, the order of calling is
  // as follows - see G4VisManager.cc.
  // DoneWithMasterThread
  // MovingToVisSubThread ) or ( SwitchToVisSubThread
  // SwitchToVisSubThread )    ( MovingToVisSubThread
  // DoneWithVisSubThread
  // MovingToMasterThread
  // SwitchToMasterThread

  // Called on the master thread before starting the vis sub-thread.
  virtual void DoneWithMasterThread ();

  // Called on the master thread after starting the vis sub-thread.
  virtual void MovingToVisSubThread ();

  // Called on the vis sub-thread at start of vis sub-thread.
  virtual void SwitchToVisSubThread ();

  // Called on the vis sub-thread when all events have been processed.
  virtual void DoneWithVisSubThread ();

  // Called on the vis sub-thread when all events have been processed.
  virtual void MovingToMasterThread ();
  
  // Called on the master thread after the vis sub-thread has terminated.
  virtual void SwitchToMasterThread ();
#endif

  //////////////////////////////////////////////////////////////
  // Access functions.
  const G4String&         GetName           () const;
  const G4String&         GetShortName      () const;
  void                    SetName           (const G4String&);
  G4int                   GetViewId         () const;
  G4VSceneHandler*        GetSceneHandler   () const;
  const G4ViewParameters& GetViewParameters        () const;
  const G4ViewParameters& GetDefaultViewParameters () const;
  G4double                GetKernelVisitElapsedTimeSeconds () const;

  virtual const std::vector<G4ModelingParameters::VisAttributesModifier>*
  GetPrivateVisAttributesModifiers() const;
  // So that privately accumulated vis attributes modifiers may be
  // concatenated with the standard vis attributes modifiers for commands
  // such as /vis/viewer/set/all and /vis/viewer/save.

  void SetViewParameters         (const G4ViewParameters& vp);
  void SetDefaultViewParameters  (const G4ViewParameters& vp);

  //////////////////////////////////////////////////////////////
  // Public utility functions.

  const G4VisAttributes*  GetApplicableVisAttributes
                                            (const G4VisAttributes*) const;

  void SetNeedKernelVisit (G4bool need);
  // Sets individual need-visit flag.

  void NeedKernelVisit ();
  // Flags all views the need to re-visit the GEANT4 kernel to refresh
  // the scene.

  void ProcessView ();
  // Used by DrawView ().  Invokes SetView ().  The basic logic is here.

protected:

  //////////////////////////////////////////////////////////////
  // Protected utility functions.

  void SetTouchable
  (const std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>& fullPath);
  // Set the touchable for /vis/touchable/set/... commands.

  void TouchableSetVisibility
  (const std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>& fullPath,
   G4bool visibility);
  // Set the touchable visibility attribute.
  // Changes the Vis Attribute Modifiers WITHOUT triggering a rebuild.

  void TouchableSetColour
  (const std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>& fullPath,
   const G4Colour&);
  // Set the touchable colour attribute.
  // Changes the Vis Attribute Modifiers WITHOUT triggering a rebuild.

  class G4Spline
  {
  public:

    // Constructors and destructor
    G4Spline();
    ~G4Spline();

    // Operations
    void AddSplinePoint(const G4Vector3D& v);
    G4Vector3D GetInterpolatedSplinePoint(G4float t);   // t = 0...1; 0=vp[0] ... 1=vp[max]
    G4int GetNumPoints();
    G4Vector3D GetPoint(int);
    // method for computing the Catmull-Rom parametric equation
    // given a time (t) and a vector quadruple (p1,p2,p3,p4).
    G4Vector3D CatmullRom_Eq(G4float t, const G4Vector3D& p1, const G4Vector3D& p2,
                             const G4Vector3D& p3, const G4Vector3D& p4);

  private:
    std::vector<G4Vector3D> vp;
    G4float delta_t;
  };
  
  //////////////////////////////////////////////////////////////
  // Data members
  G4VSceneHandler&        fSceneHandler;     // Abstract scene for this view.
  G4int            fViewId;    // Id of this instance.
  G4String         fName;
  G4String         fShortName; // Up to first ' ' character, if any.
  G4ViewParameters fVP;        // View parameters.
  G4ViewParameters fDefaultVP; // Default view parameters.
  G4double         fKernelVisitElapsedTimeSeconds = 999.;  // Default to a large number
  // Note: fKernelVisitElapsedTimeSeconds is measured in ProcessView().

  //////////////////////////////////////////////////////////////
  // Other parameters.
  G4bool           fNeedKernelVisit;  // See DrawView() for comments.
};

#include "G4VViewer.icc"

/*********************************************

Here is a minimal DrawView () as it might be implemented in the
concrete viewer.

void G4VViewer::DrawView () {  // Default - concrete view usually overrides.

  // First, a view should decide when to re-visit the G4 kernel.
  // Sometimes it might not be necessary, e.g., if the scene is stored
  // in a graphical database (e.g., OpenGL's display lists) and only
  // the viewing angle has changed.  But graphics systems without a
  // graphical database will always need to visit the G4 kernel.

  NeedKernelVisit ();  // Default is - always visit G4 kernel.
  // Note: this routine sets the fNeedKernelVisit flag of *all* the views of
  // the scene.

  ProcessView ();             // The basic logic is here.

  // Then a view may have more to do, e.g., display the graphical
  // database.  That code should come here before finally...

  FinishView ();              // Flush streams and/or swap buffers.
}

*********************************************/

#endif
