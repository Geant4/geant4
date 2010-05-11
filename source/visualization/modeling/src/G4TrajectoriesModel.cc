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
// $Id: G4TrajectoriesModel.cc,v 1.25 2010-05-11 11:21:52 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  26th August 1998.
// Model which knows how to draw GEANT4 trajectories.

#include "G4TrajectoriesModel.hh"

#include "G4ModelingParameters.hh"
#include "G4VGraphicsScene.hh"
#include "G4Event.hh"
#include "G4AttValue.hh"
#include "G4AttDef.hh"
#include "G4AttCheck.hh"

G4TrajectoriesModel::G4TrajectoriesModel ():
  fDrawingModeSet(false),
  fDrawingMode(0),
  fpCurrentTrajectory(0)
{
  fGlobalTag = "G4TrajectoriesModel for all trajectories.";
  fGlobalDescription = fGlobalTag;
}

G4TrajectoriesModel::G4TrajectoriesModel (G4int drawingMode):
  fDrawingModeSet(true),
  fDrawingMode(drawingMode),
  fpCurrentTrajectory(0) {
  fGlobalTag = "G4TrajectoriesModel for all trajectories.";
  fGlobalDescription = fGlobalTag;
}

G4TrajectoriesModel::~G4TrajectoriesModel () {}

void G4TrajectoriesModelDebugG4AttValues(const G4VTrajectory*);

void G4TrajectoriesModel::DescribeYourselfTo (G4VGraphicsScene& sceneHandler) {
  const G4Event* event = fpMP->GetEvent();
  if (event) {
    G4TrajectoryContainer* TC = event -> GetTrajectoryContainer ();
    if (TC) {
      for (G4int iT = 0; iT < TC->entries(); iT++) {
	fpCurrentTrajectory = (*TC) [iT];
	// Debug trajectory:
	// fpCurrentTrajectory->ShowTrajectory(); G4cout << G4endl;
	// Debug G4AttValues:
	// G4TrajectoriesModelDebugG4AttValues(fpCurrentTrajectory);
	if (fpCurrentTrajectory)
	  sceneHandler.AddCompound (*fpCurrentTrajectory);
      }
    }
  }
}

G4bool G4TrajectoriesModel::IsDrawingModeSet() const
{
  return fDrawingModeSet;
}

G4int G4TrajectoriesModel::GetDrawingMode() const
{
  if (fDrawingModeSet) return fDrawingMode;
  else {
    G4Exception
      ("G4TrajectoriesModel::GetDrawingMode",
       "",
       FatalException,
       "Illegal attempt to obtain i_mode."
       "\n  i_mode is an old trajectories parameter that is DEPRECATED"
       "\n  and will be removed at the next major release."
       );
    return 0;
  }
}

void G4TrajectoriesModel::SetDrawingMode(G4int drawingMode)
{
  if (fDrawingModeSet) fDrawingMode = drawingMode;
  else {
    G4Exception
      ("G4TrajectoriesModel::SetDrawingMode",
       "",
       FatalException,
       "Illegal attempt to set i_mode."
       "\n  i_mode is an old trajectories parameter that is DEPRECATED"
       "\n  and will be removed at the next major release."
       );
  }
}

// Debug material...

#include "G4VTrajectoryPoint.hh"

void G4TrajectoriesModelDebugG4AttValues(const G4VTrajectory* pTraj)
{
  std::vector<G4AttValue>* attValues = pTraj->CreateAttValues();
  if (attValues) {
    G4AttCheck attCheck(attValues, pTraj->GetAttDefs());
    G4cout << "\nProvided G4Atts:\n" << attCheck;
    if (attCheck.Check()) G4cout << "Error" << G4endl;
    else {
      std::vector<G4AttValue> standardValues;
      std::map<G4String,G4AttDef> standardDefinitions;
      attCheck.Standard(&standardValues, &standardDefinitions);
      G4cout << "\nStandard G4Atts:\n"
	     << G4AttCheck(&standardValues, &standardDefinitions);
    }
    delete attValues;
  }
  for (G4int i = 0; i < pTraj->GetPointEntries(); i++) {
    G4VTrajectoryPoint* aPoint = pTraj->GetPoint(i);
    std::vector<G4AttValue>* attValues = aPoint->CreateAttValues();
    if (attValues) {
      G4AttCheck attCheck(attValues, aPoint->GetAttDefs());
      G4cout << "\nProvided G4Atts:\n" << attCheck;
      if (attCheck.Check()) G4cout << "Error" << G4endl;
      else {
	std::vector<G4AttValue> standardValues;
	std::map<G4String,G4AttDef> standardDefinitions;
	attCheck.Standard(&standardValues, &standardDefinitions);
	G4cout << "\nStandard G4Atts:\n"
	       << G4AttCheck(&standardValues, &standardDefinitions);
      }
      delete attValues;
    }
  }	  
}
