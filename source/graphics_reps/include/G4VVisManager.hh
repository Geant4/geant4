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
// John Allison 19/Oct/1996.
// 
// Class Description:
//
// G4VVisManager is an abstract interface for the GEANT4 Visualization Manager.
// The inheritance hierarchy is:
//   G4VVisManager <- G4VisManager <- G4VisExecutive
//
// You may also write your own vis manager in place of G4VisExecutive.
//
// See example/novice/N02 to see how and when to instantiate
// G4VisExecutive (or your own vis manager).  You should *not* access
// it directly; instead you should obtain a pointer as follows:
// 
//   G4VVisManager* pVVMan =  G4VVisManager::GetConcreteInstance ();
//
// This ensures your code will link even if G4VisExecutive is not
// instantiated or even if not provided in a library.  Please protect
// your code by testing the pointer, for example, by:
//
//   if (pVVMan) pVVMan -> Draw (polyline);
//
// The Draw functions draw only "transient" objects.  This is useful
// for debugging, e.g., drawing the step in your UserSteppingAction,
// since G4Steps are not kept.
//
// Note: "permanent" objects, i.e., objects which are always
// available, such as detector geometry components, or available in an
// event after tracking has finished, such as hits, digitisations and
// trajectories, can be drawn in a transient way if you wish but it is
// usually possible to draw them in a permanent way with /vis/
// commands.  The advantage is that permanent objects can be redrawn,
// e.g., when you change view or viewer; transient objects get
// forgotten.  Also, it is possible to write a G4VUserVisAction class
// and register it to "promote" your Draw messages to "permanent" -
// see documentation.
//
// Note that the G4Transform3D argument refers to the transformation
// of the *object*, not the transformation of the coordinate syste.
//
// Note also that where a G4VisAttributes argument is specified, it
// overrides any attributes belonging to the object itself.

#ifndef G4VVISMANAGER_HH
#define G4VVISMANAGER_HH

#include "G4Transform3D.hh"
#include "G4ThreeVector.hh"       // Just a typedef Hep3Vector.
#include "G4RotationMatrix.hh"    // Just a typedef HepRotation.

class G4Polyline;
class G4Text;
class G4Circle;
class G4Square;
class G4Polymarker;
class G4Polyhedron;
class G4VSolid;
class G4VHit;
class G4VDigi;
class G4VTrajectory;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VisAttributes;

class G4VVisManager {

public: // With description

  static G4VVisManager* GetConcreteInstance ();
  // Returns pointer to actual visualization manager if a view is
  // available for drawing, else returns null.  Always check value.

public:

  G4VVisManager ();
  virtual ~G4VVisManager ();

public: // With description

  ///////////////////////////////////////////////////////////////////////
  // Draw methods for Geant4 Visualization Primitives, useful
  // for representing hits, digis, etc.

  virtual void Draw (const G4Circle&,
    const G4Transform3D& objectTransformation = G4Transform3D()) = 0;

  virtual void Draw (const G4Polyhedron&,
    const G4Transform3D& objectTransformation = G4Transform3D()) = 0;

  virtual void Draw (const G4Polyline&,
    const G4Transform3D& objectTransformation = G4Transform3D()) = 0;

  virtual void Draw (const G4Polymarker&,
    const G4Transform3D& objectTransformation = G4Transform3D()) = 0;

  virtual void Draw (const G4Square&,
    const G4Transform3D& objectTransformation = G4Transform3D()) = 0;

  virtual void Draw (const G4Text&,
    const G4Transform3D& objectTransformation = G4Transform3D()) = 0;

  ///////////////////////////////////////////////////////////////////////
  // For 2D methods, the x,y coordinates are interpreted as screen
  // coordinates, -1 < x,y < 1.  The z-coordinate is ignored.

  virtual void Draw2D (const G4Circle&,
    const G4Transform3D& objectTransformation = G4Transform3D()) = 0;

  virtual void Draw2D (const G4Polyhedron&,
    const G4Transform3D& objectTransformation = G4Transform3D()) = 0;

  virtual void Draw2D (const G4Polyline&,
    const G4Transform3D& objectTransformation = G4Transform3D()) = 0;

  virtual void Draw2D (const G4Polymarker&,
    const G4Transform3D& objectTransformation = G4Transform3D()) = 0;

  virtual void Draw2D (const G4Square&,
    const G4Transform3D& objectTransformation = G4Transform3D()) = 0;

  virtual void Draw2D (const G4Text&,
    const G4Transform3D& objectTransformation = G4Transform3D()) = 0;

  ///////////////////////////////////////////////////////////////////////
  // Draw methods for Geant4 Objects as if they were Visualization
  // Primitives.  Note that the visualization attributes needed in
  // some cases override any visualization attributes that are
  // associated with the object itself - thus you can, for example,
  // change the colour of a physical volume.

  virtual void Draw (const G4VTrajectory&) = 0;

  virtual void Draw (const G4VHit&) = 0;

  virtual void Draw (const G4VDigi&) = 0;

  virtual void Draw (const G4LogicalVolume&, const G4VisAttributes&,
    const G4Transform3D& objectTransformation = G4Transform3D()) = 0;

  virtual void Draw (const G4VPhysicalVolume&, const G4VisAttributes&,
    const G4Transform3D& objectTransformation = G4Transform3D()) = 0;

  virtual void Draw (const G4VSolid&, const G4VisAttributes&,
    const G4Transform3D& objectTransformation = G4Transform3D()) = 0;

  virtual void DrawGeometry
  (G4VPhysicalVolume*, const G4Transform3D& t = G4Transform3D());
  // Draws a geometry tree starting at the specified physical volume.

  //////////////////////////////////////////////////////////////////////
  // Optional methods that you may use to bracket a series of Draw
  // messages that have identical objectTransformation to improve
  // drawing speed.  Use Begin/EndDraw for a series of Draw messages,
  // Begin/EndDraw2D for a series of Draw2D messages.  Do not mix Draw
  // and Draw2D messages.

  virtual void BeginDraw
  (const G4Transform3D& objectTransformation = G4Transform3D()) = 0;

  virtual void EndDraw () = 0;

  virtual void BeginDraw2D
  (const G4Transform3D& objectTransformation = G4Transform3D()) = 0;

  virtual void EndDraw2D () = 0;

  //////////////////////////////////////////////////////////////////////
  // Other methods...

  virtual void GeometryHasChanged () = 0;
  // This is used by the run manager to notify a change of geometry.

  virtual void IgnoreStateChanges(G4bool);
  // This method shoud be invoked by a class that has its own event loop,
  // such as the RayTracer, material scanner, etc. If the argument is true,
  // the following state changes among Idle, GeomClosed and EventProc are
  // caused by such a class, and thus not by the ordinary event simulation.
  // The same method with false should be invoked once such an event loop
  // is over.

  virtual void NotifyHandlers () {}
  // Notify scene handlers (G4VGraphicsScene objects) that the scene
  // has changed so that they may rebuild their graphics database, if
  // any, and redraw all views.

  virtual void DispatchToModel(const G4VTrajectory&) = 0;
  // Draw the trajectory.

  virtual G4bool FilterTrajectory(const G4VTrajectory&) = 0;
  // Trajectory filter

  virtual G4bool FilterHit(const G4VHit&) = 0;
  // Hit filter

  virtual G4bool FilterDigi(const G4VDigi&) = 0;
  // Digi filter

#ifdef G4MULTITHREADED

  virtual void SetUpForAThread() = 0;
  // This method is invoked by G4WorkerRunManager

#endif

protected:

  static void SetConcreteInstance (G4VVisManager*);

  static G4VVisManager* fpConcreteInstance;  // Pointer to real G4VisManager.

};

#endif
