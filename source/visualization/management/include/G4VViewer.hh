// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VViewer.hh,v 1.6 2000-02-21 16:05:21 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

class G4VSceneHandler;

class G4VViewer {

public: // With description

  friend G4std::ostream& operator << (G4std::ostream& os, const G4VViewer& v);

  G4VViewer (G4VSceneHandler& scene, G4int id, const G4String& name = "");
  virtual ~G4VViewer ();

  // For G4RWTPtrOrderedVector...
  G4bool operator == (const G4VViewer& view) const;

  //////////////////////////////////////////////////////////////
  // View manipulation functions.

  virtual void SetView () = 0;
  // Take view parameters and work out model/view transformation,
  // projection transformation, lighting, etc.

  virtual void ClearView () = 0;
  // Clear screen/viewing buffers.

  virtual void DrawView () = 0;
  // Draw scene - see example of a minimal function at end of this file.

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

  //////////////////////////////////////////////////////////////
  // Access functions.
  const G4String&         GetName           () const;
  void                    SetName           (const G4String&);
  G4int                   GetViewId         () const;
        G4VSceneHandler*         GetScene          () const;
  const G4ViewParameters& GetViewParameters () const;
  void SetViewParameters  (const G4ViewParameters& vp);

  //////////////////////////////////////////////////////////////
  // Public utility functions.

  const G4VisAttributes*  GetApplicableVisAttributes
                                            (const G4VisAttributes*) const;

  void SetNeedKernelVisit ();
  // Sets individual need-visit flags.

  void NeedKernelVisit ();
  // Flags all views the need to re-visit the GEANT4 kernel to refresh
  // the scene.

protected:

  void ProcessView ();
  // Used by DrawView ().  Invokes SetView ().  The basic logic is here.

  //////////////////////////////////////////////////////////////
  // Data members
  G4VSceneHandler&        fSceneHandler;     // Abstract scene for this view.
  G4int            fViewId;    // Id of this instance.
  G4String         fName;
  G4String         fShortName; // Up to first ' ' character, if any.
  G4ViewParameters fVP;        // Viewing parameters.

  //////////////////////////////////////////////////////////////
  // Other parameters.
  G4bool           fModified;         // If View Parameters have been modified.
  G4bool           fNeedKernelVisit;  // See DrawView() for comments.
};

#include "G4VViewer.icc"

/*********************************************

Here is a minimal DrawView () as it might be implemented in the
concrete view.

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
