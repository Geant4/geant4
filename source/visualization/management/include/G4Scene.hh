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
// $Id: G4Scene.hh,v 1.17 2006/06/29 21:27:56 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// 
// Scene  John Allison  19th July 1996.
//
// Class Description:
//
// Defines the viewable scene.

#ifndef G4SCENE_HH
#define G4SCENE_HH

#include "globals.hh"
#include "G4ios.hh"

class G4VPhysicalVolume;

#include "G4VisExtent.hh"
#include "G4Point3D.hh"
#include "G4VModel.hh"
#include <vector>

class G4Scene {

public: // With description

  friend std::ostream& operator << (std::ostream& os, const G4Scene& d);

  enum {UNLIMITED = -1};

  G4Scene (const G4String& name = "scene-with-unspecified-name");
  ~G4Scene ();

  // Makes use of default (compiler generated) copy constructor and
  // assignment operator.

  G4bool operator == (const G4Scene&) const;
  G4bool operator != (const G4Scene&) const;

  //////////////////////////////////////////////////////
  // Get functions...

  const G4String& GetName () const;

  G4bool IsEmpty () const;

  const std::vector<G4VModel*>& GetRunDurationModelList () const;
  // Contains models which are expected to last for the duration of
  // the run, for example geometry volumes.

  const std::vector<G4VModel*>& GetEndOfEventModelList () const;
  // Contains models which are described at the end of event when the
  // scene is current.

  const G4VisExtent& GetExtent () const;
  // Overall extent of all objects in the scene.

  const G4Point3D& GetStandardTargetPoint () const;
  // Usually centre of extent.  See G4ViewParameters for definition.

  G4bool GetRefreshAtEndOfEvent () const;
  // If true, the visualization manager will request viewer to refresh
  // "transient" objects, such as hits, at end of event.  Otherwise
  // they will be accumulated.

  G4bool GetRefreshAtEndOfRun () const;
  // If true, the visualization manager will request viewer to refresh
  // "transient" objects, such as hits, at end of run.  Otherwise
  // they will be accumulated.

  G4bool GetRecomputeTransients () const;
  // If true, the visualization manager will attempt to re-run events
  // in order to reconstruct trajectories, etc., when the viewer
  // requests a scene traversal.

  //////////////////////////////////////////////
  // Add and Set functions...

  G4bool AddRunDurationModel (G4VModel*, G4bool warn = false);
  // Adds models of type which are expected to last for the duration
  // of the run, for example geometry volumes.
  // Returns false if model is already in the list.
  // Prints warnings if warn is true.

  G4bool AddWorldIfEmpty (G4bool warn = false);
  // In some situations, if the user asks for a drawing and has not
  // yet set any run duration models it makes sense to put the "world"
  // in there by default.
  // Returns false if model is already in the list.
  // Prints warnings if warn is true.

  G4bool AddEndOfEventModel (G4VModel*, G4bool warn = false);
  // Adds models of type which are described at the end of event when
  // the scene is current.
  // Returns false if model is already in the list.
  // Prints warnings if warn is true.

  std::vector<G4VModel*>& SetRunDurationModelList ();
  // Allows you to change the model list - do with care!

  std::vector<G4VModel*>& SetEndOfEventModelList ();
  // Allows you to change the model list - do with care!

  void SetRefreshAtEndOfEvent(G4bool);
  // If set true, the visualization manager will request viewer to
  // refresh "transient" objects, such as hits, at end of event.
  // Otherwise they will be accumulated.

  void SetRefreshAtEndOfRun(G4bool);
  // If set true, the visualization manager will request viewer to
  // refresh "transient" objects, such as hits, at end of run.
  // Otherwise they will be accumulated.

  void SetRecomputeTransients(G4bool);
  // If set true, the visualization manager will attempt to re-run
  // events in order to reconstruct trajectories, etc., when the
  // viewer requests a scene traversal.

  //////////////////////////////////////////////
  // Other functions...

  void CalculateExtent();
  // (Re-)calculates the extent from the extents of its models.

private:
  G4String fName;
  std::vector<G4VModel*> fRunDurationModelList;
  std::vector<G4VModel*> fEndOfEventModelList;
  G4VisExtent fExtent;
  G4Point3D   fStandardTargetPoint;
  G4bool      fRefreshAtEndOfEvent;
  G4bool      fRefreshAtEndOfRun;
  G4bool      fRecomputeTransients;
};

#include "G4Scene.icc"

#endif
